%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Math Lead & Secondary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Provided by the Herzenberg Lab at Stanford University 
%   License: BSD 3 clause
%
classdef WebDownload
    properties(Constant)
        %change this to true for white box testing of google drive
        SIMULATE_GOOGLE_DRIVE_FALL_BACK=false;
        
        %change this if the fall back google drive is not stephen meehan's
        DEFAULT_GOOGLE_DRIVE_DIRECTORY='https://drive.google.com/file/d/1r_ByRnrq7SWvfLa2uGdwE31IY1ClXVdy/view?usp=sharing';
        
        
        PROPERTY_GOOGLE_DRIVE_DIRECTORY='googleDriveDirectory';
        
        HOST='http://cgworkspace.cytogenie.org/';
        PATH='GetDown2/demo';
        DOCUMENTS_FOLDER_FOR_EXAMPLES='run_umap/examples';
        
        KEY_UMAP_EXAMPLES='umapExamples';
        KEY_DEMO='demo';
        KEY_AutoGate='facs';
        
        PROPERTY_HOST='host';
        DEFAULT_HOST='http://cgworkspace.cytogenie.org';
        DEFAULT_BAD_HOST='http://cgworkspace.cytogenies.org';
        HOST_FILE='.herzenbergLabHosts.properties';        
    end
    
    methods(Static)
        
        function UpdateGoogleDirProperties(fileSpecOnGoogleDrive, fileName)
            folder=fullfile(File.Home, 'Google Drive');
            spec=qualify(fileSpecOnGoogleDrive, 'file spec');
            if nargout<2
                fileName=fullfile(folder, 'dir.properties');
            end
            fileName=qualify(fileName, 'directory property file');
            m=File.ReadProperties(fileName);
            subFldr=fileparts(spec(length(folder)+2:end));
            fs=dir(spec);
            gd='^(?<size>\d+)\s+(?<url>https://drive.google.com/.*)';
            gd2='^(?<size>\d+)\s+';
            N=length(fs);
            if N==0
                whine=Html.Wrap(['File does not exist<br>'...
                    Html.WrapBoldSmall(String.ToHtml(spec))]);
                msg(whine);
                warning(whine);
                return;
            end
            WebDownload.StateSharingRequirements
            
            changes=0;
            for i=1:N
                f=fs(i);
                if ~startsWith(f.name, '.') && f.bytes>0
                    key=fullfile(subFldr, f.name);
                    if ispc
                        key=strrep(key, filesep, '/');
                    end
                    value=num2str(f.bytes);
                    if ~m.containsKey(key)
                        disp(['      Adding to Google dir "' key '"=' value]);
                        m.put(key, [value ' ']);
                        changes=changes+1;
                    else
                        old=m.get(key);
                        flds=regexp(old, gd, 'names');
                        if isempty(flds)
                            flds=regexp(old, gd2, 'names');
                            if ~isempty(flds)
                                flds.url='';
                            end
                        end
                        if isempty(flds)
                            warning(['Bad/incomplete format ' key '=' old ]);
                            disp(['      Correcting Google dir "' key '"=' value]);
                            m.put(key, [value ' ']);
                            changes=changes+1;
                        else
                            if ~strcmp(flds.size, value)
                                disp([' .     Changing Google dir size from ' ...
                                    flds.size ' to ' value ' for "' key '"']);
                                m.put(key, [value ' ' flds.url]);
                                changes=changes+1;
                            else
                                disp(['       Google dir unchanged for "'...
                                    key '"="' old '"']);
                            end
                        end
                    end
                end
            end
            if changes>0
                File.SaveProperties(m, fileName, true);
            end
            function qualified=qualify(fileThing, word)
                jf=java.io.File(fileThing);
                
                if ~jf.isAbsolute
                    qualified=fullfile(folder, fileThing);
                else
                    qualified=fileThing;
                end
                assert(startsWith(qualified, folder), ...
                    ['The Google ' word ...
                    ' must be on your Google Drive!']);
            end
        end
        
        function [hosts, googleDir, m]=GetHosts(urlPath)
            testingGoogleDrive=WebDownload.SIMULATE_GOOGLE_DRIVE_FALL_BACK;
            fileName=fullfile(File.Home, WebDownload.HOST_FILE);
            if ~testingGoogleDrive
               	m=JavaProperties(fileName);
            else
                m=JavaProperties;
            end
            if ~testingGoogleDrive
                m0=m.addIfMissing(WebDownload.PROPERTY_HOST, WebDownload.DEFAULT_HOST);
            else
                m0=m.addIfMissing(WebDownload.PROPERTY_HOST, WebDownload.DEFAULT_BAD_HOST);
            end
            m1=m.addIfMissing(WebDownload.PROPERTY_HOST, 'https://drive.google.com');
            m2=m.addIfMissing(WebDownload.PROPERTY_HOST, 'http://54.39.2.45');
            m3=m.putIfMissing(WebDownload.PROPERTY_GOOGLE_DRIVE_DIRECTORY, ...
                WebDownload.DEFAULT_GOOGLE_DRIVE_DIRECTORY);
            if m0 || m1 || m2 || m3
                if ~testingGoogleDrive
                    m.save(fileName);
                end
            end
            specificKey=[WebDownload.PROPERTY_HOST '.' urlPath];
            hosts0=m.getAll(specificKey);
            specificKey=[WebDownload.PROPERTY_GOOGLE_DRIVE_DIRECTORY '.' urlPath];
            googleDir=m.get(specificKey);
            hosts=[hosts0 m.getAll(WebDownload.PROPERTY_HOST)];
            if isempty(googleDir)
                googleDir=m.get(WebDownload.PROPERTY_GOOGLE_DRIVE_DIRECTORY);
            end
        end
        
        function StateSharingRequirements
            disp('Any url like https://drive.google.com/run_umap/examples/sample10k.csv');
            disp('  gets converted to an internal Google Drive export format  ');
            disp('  following the logic at https://www.wonderplugin.com/online-tools/google-drive-direct-link-generator/');
            disp('  This conversion requires a properties file ');
            disp('  that records the relative filepath and size and the sharing');
            disp('  link produced by the Google Drive user interface');
            disp('  ');
            disp('THEN code like WebDownload.java can directly download the ');
            disp('  URL in export format ONLY if both the file ');
            disp('  and the parent folder are shared read-only to everyone!!');
            disp(' ');
            disp(' This does not work for files like tar that are too big to download');
        end
        
        function [url, googleDir]=ResolveHost(urlPath)
            app=BasicMap.Global;
            map=app.urlMap;
            [hosts, googleDir]=WebDownload.GetHosts(urlPath);
            N=length(hosts);
            for i=1:N
                beforePath=hosts{i};
                [host,path,port]=WebDownload.UrlParts(beforePath);
                if isempty(port) || port<=0 
                    port=80;
                end
                badKey=['bad:' host ':' num2str(port)];
                if ~map.containsKey(badKey)
                    [ok, issues]=WebDownload.CanReach(host, port);
                    if ok
                        url=WebDownload.FullUrl(beforePath, urlPath);
                        return;
                    end
                    map.put(badKey, badKey);
                end
            end
            if N>0
                map.clear;
                url=urlPath;
            else
                url='';
            end
        end
        
        function [ok, issues]=CanReach(host, port, timeout)
            if nargin<3
                timeout=4000;
            end
            if isempty(port) || port<=0 
                port=80;
            end
            issues=java.util.ArrayList;
            ok=edu.stanford.facs.swing.WebDownload.CanReachHost(host, port, timeout, issues);
        end
        
        function [url, googleDir]=ResolveUrl(file, key)
            if nargin<2
                key='run_umap/examples';
                if nargin<1
                    file='';
                end
            end
            app=BasicMap.Global;
            map=app.urlMap;
            firstTime=map.size==0;
            url=[];
            googleDir=[];
            if WebDownload.SIMULATE_GOOGLE_DRIVE_FALL_BACK
                map.clear;
            end
            if ~map.containsKey(key)
                [url, googleDir]=WebDownload.ResolveHost(key);
                map.put(key, url);
            else
                url=map.get(key);
            end
            if ~firstTime
                [host,path,port]=WebDownload.UrlParts(url);
                if ~WebDownload.CanReach(host, port)
                    map.clear;
                    [url, googleDir]=WebDownload.ResolveUrl(file, key);
                end
            end
            if ~isempty(url)
                url=[url '/' strrep(file, '\', '/')];
            end
        end
        
        function url=Url(file, path, host)
            if nargin<3
                host=WebDownload.HOST;
                if nargin<2
                    path=WebDownload.PATH;
                    if nargin<1
                        file='';
                    end
                end
            end
            url=[host path '/' file];
        end
        
        function [fullFile, exists, downloaded]=Download(...
                fromUrl, whereOnScreen, toLocalFolder, allowCancel)
            if nargin<4
                allowCancel=true;
                if nargin<3
                    toLocalFolder=tempdir;
                    if nargin<2
                        whereOnScreen='center';         
                    end
                end
            end
            exists=false;
            downloaded=false;
            [host,path]=WebDownload.UrlParts(fromUrl);
            if isempty(host)
                fullFile=[];
                msgError(Html.WrapHr(['Not a URL...<br>'...
                    Html.WrapSmallTags(['(Bad format "' ...
                    fromUrl '")' ])]));
                return;
            end
            if isempty(path)
                remoteFile='tempWebDownloadFile';
            else
                [~, f, e]=fileparts(path);
                remoteFile=[f e];
            end
            fullFile=fullfile(toLocalFolder, remoteFile);
            if exist(fullFile, 'file') && ~isequal(tempdir, toLocalFolder)
                if ~askYesOrNo(Html.WrapHr(['Overwrite prior ?<br>' ...
                        Html.WrapSmallTags(['(' fullFile ')'])]), ...
                        'Exists  ....', whereOnScreen, true)
                    exists=true;
                    return;
                end
            end
            [cancelledByUser, bad]=WebDownload.Get(...
                {fromUrl}, {fullFile}, false, allowCancel, whereOnScreen); 
            downloaded=~cancelledByUser && ~bad;
            fprintf('Downloaded=%d, cancelled=%d and bad=%d\n', ...
                downloaded, cancelledByUser, bad);
            if ~downloaded
                fullFile=[];
            end                
        end
        
        function [bytesPerUrl, didUrlConnect, problemsPerUrl, dwl]=GetSize(urls)
            try
                dwl=javaObjectEDT('edu.stanford.facs.swing.WebDownload');
            catch ex
                jar=fullfile(fileparts(mfilename('fullpath')), 'webDownload.jar');
                javaaddpath(jar);
                try
                    dwl=javaObjectEDT('edu.stanford.facs.swing.WebDownload');
                catch ex
                    ex.getReport
                    error('Trouble using WebDownload');
                end
            end
            dwl.waitWhenDone=false;
            dwl.allowCancel=false;
            N=length(urls);
            bytesPerUrl=zeros(1,N);
            didUrlConnect=false(1,N);
            problemsPerUrl=cell(1,N);
            for i=1:N
                problemsPerUrl{i}=java.util.ArrayList;
            	con=dwl.connect(urls{i}, problemsPerUrl{i});
                if ~isempty(con.http)
                    didUrlConnect(i)=true;
                    bytesPerUrl(i)= con.size;
                    con.disconnect();
                end
						
            end
        end
        
        
        %where only used if Gui.m is present
        function [cancelledByUser, bad, dwl, dlg]=...
                Get(urls, localFiles,  waitWhenDone, allowCancel, where)
            if nargin<5
                where='center';
                if nargin<4
                    allowCancel=true;
                    if nargin<3
                        waitWhenDone=true;
                    end
                end
            end
            try
                dwl=javaObjectEDT('edu.stanford.facs.swing.WebDownload');
            catch ex
                jar=fullfile(fileparts(mfilename('fullpath')), 'webDownload.jar');
                javaaddpath(jar);
                try
                    dwl=javaObjectEDT('edu.stanford.facs.swing.WebDownload');
                catch ex
                    ex.getReport
                    error('Trouble using WebDownload');
                end
            end
            dlg=javaObjectEDT(dwl.dlg);
            try
                Gui.LocateJava(dlg, ...
                    Gui.JWindow(get(0, 'CurrentFigure')), where);
            catch
                WebDownload.LocateJavaOnScreen(dlg, where);
            end
            dwl.waitWhenDone=waitWhenDone;
            dwl.allowCancel=allowCancel;
            cancelledByUser=~dwl.go(urls,localFiles, ~isempty(where));
            bad=dwl.bad;
        end
        
        function txt=ReadText(url)
            localFile=[tempname '.txt'];
            [cancelledByUser, bad]=...
                WebDownload.Get({url}, {localFile},  false, false, []);
            if ~cancelledByUser && ~bad
                txt=File.ReadTextFile(localFile);
            else
                txt=[];
            end
        end
        %where only used if Gui.m is present
        function [ok, cancelledByUser]=GetAndUnzip(url, zipFile, waitWhenDone, ...
                allowCancel, where)
            if nargin<5
                where='center';
                if nargin<4
                    allowCancel=true;
                    if nargin<3
                        waitWhenDone=true;
                    end
                end
            end
            ok=false;
            [cancelledByUser, bad, dwl, dlg]=WebDownload.Get({url}, ...
                {zipFile}, waitWhenDone, allowCancel, where);
            if ~cancelledByUser && bad==0
                dlg.setModal(false);
                dwl.progressBar.setValue(0);
                dlg.setVisible(true);
                dwl.south.setText('Unzipping file now');
                try
                    zipFldr=fileparts(zipFile);
                    if isempty(zipFldr)
                        unzip(zipFile);
                    else
                        unzip(zipFile, zipFldr);
                    end
                    dwl.progressBar.setValue(dwl.progressBar.getMaximum);
                    MatBasics.RunLater(@(h,e)quiet,2);
                    ok=true;
                catch ex
                    ex.getReport
                end
            end
            delete(zipFile);
                    
            function quiet
                dlg.setVisible(false);
            end
        end
        
        function [x,y]=LocateJavaOnScreen(javaComponent, where)
            size=javaComponent.getSize;
            if size.width==0 || size.height==0
                size=javaComponent.getPreferredSize;
            end
            if nargin<2
                where='center';
            end
            p=get(0, 'ScreenSize');
            where=strrep(where,'-', '');
            where=strrep(where,'+','');
            [x, y]=WebDownload.LocateWidthHeight(true, ...
                size.width, size.height, p(1), ...
                p(2), p(3), p(4), where);
            javaComponent.setLocation(x, y);
        end
        
        function [x,y]=LocateWidthHeight(isTop0, width, height, ...
                refX, refY, refWidth, refHeight, where)
            if isempty(where)
                where='center';
            end
            w=String(lower(where));            
            hCenter=~w.contains('east') && ~w.contains('west');
            vCenter=~w.contains('north') && ~w.contains('south');
            centerX=refX+refWidth/2-width/2;
            centerY=refY+refHeight/2-height/2;
            east=w.contains('east');
            top=w.contains('north');
            plusPlus=w.contains('++');
            plus=~plusPlus&&w.contains('+');
            if hCenter && vCenter
                x=centerX;
                y=centerY;
            else
                if width>refWidth
                    if east
                        W=refWidth*.75;
                    else
                        W=refWidth*1.3;
                    end
                else
                    W=width;
                end
                if east
                    if plus
                        h=(refX+refWidth)-W*.4;
                    elseif plusPlus
                        h=(refX+refWidth)-W*.1;
                    else
                        h=(refX+refWidth)-W;
                    end
                else
                    if plus
                        h=refX-(W*.6);
                    elseif plusPlus
                        h=refX-(W*.9);
                    else
                        h=refX;
                    end
                end
                if (top && isTop0) || (~top && ~isTop0)
                    if plus
                        v=refY-height*.25;
                    elseif plusPlus
                        v=refY-height*.9;
                    else
                        v=refY;
                    end
                elseif (top && ~isTop0) || (~top && isTop0)
                    if plus
                        v=(refY+refHeight)-height*.66;
                    elseif plusPlus
                        v=(refY+refHeight)-height*.01;                    
                    else
                        v=(refY+refHeight)-height;
                    end
                end
                if vCenter
                    x=h;
                    y=centerY;
                elseif hCenter
                    x=centerX;
                    y=v;
                else
                    x=h;
                    y=v;
                end
            end
        end
        
        function file=GetZipIfMissing(file, url)
            if ischar(file) && ~exist(file, 'file')
                [fldr,fn]=fileparts(file);
                zipFileName=[fn '.zip'];
                if nargin<2
                    url=fullfile(WebDownload.HOST, WebDownload.PATH);
                end
                url=WebDownload.FullUrl(url, zipFileName);
                zipFile=fullfile(fldr, zipFileName);
                [ok, cancelledByUser]=WebDownload.GetAndUnzip(url, ...
                    zipFile, false, true, 'center');
                if cancelledByUser
                     msg(Html.WrapHr(['Cancelling leaves the required '...
                         'file missing...<br>' BasicMap.Global.smallStart...
                         '<b>' file '</b>' BasicMap.Global.smallEnd]));
                     file=[];
                elseif ~ok
                    msg(Html.WrapHr(['Required file is missing...<br>' ...
                        BasicMap.Global.smallStart...
                        '<b>' file '</b>' BasicMap.Global.smallEnd]));
                    file=[];
                end
            end
        end
        
        function url=FullUrl(beforePath, afterPath)
            if endsWith(beforePath, '/')
                beforePath=beforePath(1:end-1);
            end
            if startsWith(afterPath, '/')
                afterPath=afterPath(2:end);
            end
            url=[beforePath '/' afterPath];
        end
        
        function [host, path, port, query, protocol]=UrlParts(url)
            u=matlab.net.URI(url);
            path=char(u.EncodedPath);
            host=char(u.Host);
            query=char(u.EncodedQuery);
            idx=String.IndexOf(url, ':');
            if idx>0
                protocol=url(1:idx-1);
            else
                protocol='';
            end
            port=u.Port;
        end

        function [fullFile, exists, downloaded, isDifferent]...
                =GetFile(fileName, localFolder, remoteFolder, ...
                force, complainEvenIfExists)
            if nargin<5
                complainEvenIfExists=true;
                if nargin<4
                    force=false;
                end
            end
            fullFile=fullfile(localFolder, fileName);
            File.mkDir(localFolder);
            downloaded=false;
            isDifferent=false;
            exists=exist(fullFile, 'file');
            if force || ~exists
                if nargout>3
                    bak=tempname;
                    copyfile(fullFile, bak);
                end
                url=WebDownload.ResolveUrl(fileName, remoteFolder);
                [cancelled, bad]=WebDownload.Get({url}, {fullFile}, false, false, 'south');
                downloaded= ~cancelled && ~bad;
                if ~downloaded
                    if ~exists || complainEvenIfExists
                        msg(['<html>Could not download the file "' ...
                            fileName '"<br>'...
                            Html.WrapColor('From:  ', 'blue')...
                            Html.WrapSmallTags(url) '<br>' ...
                            Html.WrapColor('To:  ', 'blue') ...
                            Html.WrapSmallTags(fileparts(fullFile)) '<hr></html>']);
                    end
                else
                    exists=true;
                    if nargout>3
                        isDifferent=File.Diff(bak, fullFile);
                    end
                end
            end
        end
        
        function fldr=LocalExamplesFolder
            fldr=fullfile(File.Home, 'Documents', ...
                WebDownload.DOCUMENTS_FOLDER_FOR_EXAMPLES);
            File.mkDir(fldr);
        end
        
        function file=GetExampleIfMissing(file, localExamplesFolder)
            if ~exist(file, 'file')
                p=fileparts(file);
                if isempty(p)
                    if nargin<2
                        localExamplesFolder=WebDownload.LocalExamplesFolder;
                    end
                    file=fullfile(localExamplesFolder, file);
                    if ~exist(file, 'file')
                        csvFiles=WebDownload.RelocateExamples(...
                            {file}, true, {}, localExamplesFolder);
                        file=csvFiles{1};
                    end
                end
            end
        end
        
        function [result, existence, missingFiles]=...
                RelocateExamples(files, tryDownload, ignore, localExamplesFolder)
            if nargin<4
                localExamplesFolder=WebDownload.LocalExamplesFolder;
                if nargin<3
                    ignore={};
                    if nargin<2
                        tryDownload=true;
                    end
                end
            end
            testExistence=nargout>1;
            existence=[];
            missingFiles={};
            argType='cell';
            if isstruct(files)
                argType='struct';
                args=files;
                files={};
                if ischar(args.csv_file_or_data)
                    files{end+1}=args.csv_file_or_data;
                end
                if ~isempty(args.label_file)
                    files{end+1}=args.label_file;
                end
                if ~isempty(args.template_file)
                    files{end+1}=args.template_file;
                end
                if ~isempty(args.color_file)
                    files{end+1}=args.color_file;
                end
            elseif ischar(files)
                argType='char';
                files={files};
            end
            missingExamples=java.util.HashMap;
            N=length(files);
            fileUrls={};
            localFiles={};
            for i=1:N
                file=files{i};
                if ~exist(file, 'file')
                    if ~ispc
                        if file(1)=='~'
                            file=[File.Home file(2:end)];
                        end
                    end
                end
                if ~exist(file, 'file')
                    [fldr, fn, ext]=fileparts(file);
                    if isempty(fldr) || isequal(localExamplesFolder, fldr)
                        fldr=localExamplesFolder;
                        File.mkDir(fldr);
                        file2=fullfile(fldr, [fn ext]);
                        missingExamples.put(java.lang.String(file), java.lang.String(file2));
                        
                        if ~exist(file2, 'file') && tryDownload
                            fileUrls{end+1}=...
                                WebDownload.ResolveUrl([fn ext]);
                            localFiles{end+1}=file2;
                        end
                    end
                end
            end
            downloadFailures=[];
            if ~isempty(fileUrls)
                nMissing=length(fileUrls);
                [cancelledByUser, bad]=WebDownload.Get(...
                    fileUrls, localFiles, false);
                if cancelledByUser
                elseif bad==0
                    msg(Html.WrapHr([ String.Pluralize2('file', ...
                        nMissing) ' downloaded to<br><b>' ...
                        BasicMap.Global.smallStart ...
                        localExamplesFolder...
                        BasicMap.Global.smallEnd '</b>']), 8, ...
                        'south east+');
                else
                    
                    if bad==nMissing
                        downloadFailures=[String.Pluralize2('file', ...
                            nMissing) ' NOT downloaded'];
                    else
                        downloadFailures=[ num2str(bad) ' of '...
                            String.Pluralize2('file', ...
                            nMissing) ' NOT downloaded'];
                    end
                end
            end
            if isequal(argType, 'struct')
                if ischar(args.csv_file_or_data)
                    args.csv_file_or_data=grab(args.csv_file_or_data);
                end
                args.label_file=grab(args.label_file);
                args.template_file=grab(args.template_file);
                args.color_file=grab(args.color_file);
            else
                for i=1:N
                    files{i}=grab(files{i});
                end
            end
            if isequal(argType, 'cell')
                result=files;
            elseif isequal(argType, 'struct')
                result=args;
            else %argType == char
                result=files{1};
            end
            if any(~existence) && ~isempty(missingFiles)
                if ~isempty(ignore)
                    [~, canIgnore]=StringArray.EndsWith(missingFiles, ignore);
                else
                    canIgnore=false;
                end
                if ~canIgnore
                    app=BasicMap.Global;
                    html=Html.Wrap([app.h2Start 'Missing files' ...
                        app.h2End downloadFailures app.smallStart ...
                        Html.ToList(missingFiles, 'ol') app.smallEnd '<hr>']);
                    msgWarning(html, 11, 'south', 'Missing files...');
                else
                    existence(:)=2;
                end
            end
            
            function out=grab(in)
                out=in;
                if ~isempty(in)
                    key=java.lang.String(in);
                    if missingExamples.containsKey(key)
                        if tryDownload
                            out=char(missingExamples.get(key));
                        else
                            out='';
                        end
                    else
                        if in(1)=='~'
                            out=[File.Home in(2:end)];
                        end
                    end
                end
                if testExistence
                    existence(end+1)=exist(out, 'file');
                    if ~existence(end)
                        missingFiles{end+1}=out;
                    end
                end
            end
        end
    end
end