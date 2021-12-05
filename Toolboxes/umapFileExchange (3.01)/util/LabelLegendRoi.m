classdef LabelLegendRoi < handle
    properties(SetAccess=private)
        xyData;
        javaLegend;
        btnCreate;
        btns;
        btnLbls;
        map;
        asked=false;
        buildingMultiplePolygons=false;
        dfltFile;
        save_output=false;
        output_folder=[];
        ax;
        roiCallback
    end
    
    properties
        labels;    
        percentClosest=.94;
    end
    
    methods
        function this=LabelLegendRoi(xyData, roiCallback, ...
                javaLegend, btns, btnLbls, labels)
            this.xyData=xyData;
            this.roiCallback=roiCallback;
            this.javaLegend=javaLegend;
            this.btns=btns;
            this.btnLbls=btnLbls;
            this.map=Map;
            if nargin>5
                this.labels=labels;
            end
        end
        
        function followParentFig(this, parentFig, where, close)
            SuhWindow.Follow(this.javaLegend, parentFig, where, close);
        end
        
        function [lbls, btns]=getSelected(this)
            [lbls, btns]=QfTable.GetSelectedLabels(this.btns, this.btnLbls);
        end
        
        
        function processSelections(this, grouped)
            assert(size(this.xyData, 1)==length(this.labels));
            [lbls, btns_]=this.getSelected;
            if grouped
                this.makeOrDeleteOne(lbls);
            else
                nChoices=length(lbls);
                for i=1:nChoices
                    this.buildingMultiplePolygons=i<nChoices;
                    this.makeOrDeleteOne(lbls(i), btns_{i});
                    drawnow;
                end
            end
        end
        
        function roi=makeOrDeleteOne(this, lbls, btn)
            key=num2str(sort(lbls));
            if this.map.containsKey(key)
                roi=this.map.remove(key);
                if ~isempty(roi)
                    delete(roi);
                    this.map.remove([key '.name']);
                end
                return;
            end
            if length(lbls)>1
                idxs=ismember(this.labels, lbls);
            else
                idxs=this.labels==lbls;
            end
            clump=this.xyData(idxs,:);
            [~,idxs]=MatBasics.PercentCloses(clump, this.percentClosest);
            [roi, strXy]=RoiUtil.NewForXy(...
                this.ax, clump(idxs,:), .00121, this.roiCallback);
            this.map.set(key, roi);
            if nargin>2
                N=length(lbls);
                if N>1
                    [name, cancelled]=inputDlg(...
                        'Enter summary name for properties', ...
                        'ROI name', key);
                    if cancelled
                        name=key;
                    end
                else
                    name=LabelLegendRoi.Name(btn);
                end
                this.map.set([key '.name'], name);
            end
        end
        
        function save(this)
            if LabelLegendRoi.Save(this.map, this.save_output,...
                    this.output_folder, this.dfltFile, this.asked)
                this.asked=true;
            end
        end
        
        function updateLegendGui(this, ax, showCreateBtn, parentFig, where, close)
            this.ax=ax;
            if showCreateBtn
                jcmp=this.javaLegend.getContentPane;
                jp=Gui.Panel;
                btn=Gui.NewBtn(...
                    Html.WrapSmallBold('Regions of interest'), ...
                    @(h,e)dropDownMenu(this, h), ...
                    'Create/remove polygons for above subset selection', ...
                    'polygonGate.png');
                jp.add(btn);
                jcmp.add(jp, 'South')
                this.javaLegend.pack;
            end
            drawnow;
            this.followParentFig(parentFig, where, close)
        end
       
        function setSaveInfo(this, dfltFile, save_output, output_folder)
            this.dfltFile=dfltFile;
            this.save_output=save_output;
            this.output_folder=output_folder;
        end
        
    end
    
    methods(Access=private)
        
        function dropDownMenu(this, h)
            jm=PopUp.Menu;
            app_=BasicMap.Global;
            [lbls, btns_]=this.getSelected;
            nSelectedLbls=length(lbls);
            c={};
            enable=false;
            c{1}=Gui.NewMenuLabel(jm,'ROI options (regions of interest)');
            jm.addSeparator;
            callbackForIndividualOnes=...
                @(h,e)processSelections(this, false);
            if nSelectedLbls>1
                c{2}=Gui.NewMenuItem(jm, sprintf(['Create/remove '...
                    'individual polygons for above %d selections'],...
                    nSelectedLbls), callbackForIndividualOnes);
                enable=true;
            else
                if nSelectedLbls==0
                    word='NOTHING selected!!';
                else
                    word=Html.StripHtmlWord(char(btns_{1}.getText));
                end
                c{2}=Gui.NewMenuItem(jm,...
                    sprintf('<html>Create/remove polygon for "%s"</html>', ...
                    word), callbackForIndividualOnes, [],[],...
                    nSelectedLbls==1);
            end
            c{end+1}=Gui.NewMenuItem(jm, sprintf(...
                'Create/remove one single polygon for above %d selections', nSelectedLbls), ...
                @(h,e)processSelections(this, true),[],[],enable);
            c{end+1}=Gui.NewMenuItem(jm, ...
                'Save polygon information', @(h,e)saveLegendRois());
            jm.show(h, 15, 15);
            
            function saveLegendRois
                if this.map.size<1
                    msg('Nothing to save ...');
                else
                    this.save;
                end
            end
        end
    end
    
    methods(Static)
        function name=Name(btn)
            name=char(edu.stanford.facs.swing.Basics.RemoveXml(btn.getText));
            name=strrep(name, '&bull;', '');
        end
        
        function ok=Save(map, save_output, output_folder, dfltFile, asked)
            ok=false;
            
            if map.size<1
                msg('No named regions to save yet...');
                return;
            end
            if ~save_output
                [folder, file]=uiPutFile(output_folder, dfltFile, ...
                    BasicMap.Global, 'LabelLegendRoi.save', ...
                    'Save labeled ROI properties');
                if isempty(folder)
                    return;
                end
                file=fullfile(folder,file);
            else
                file=fullfile(output_folder, dfltFile);
            end
            props=java.util.Properties;
            keys=map.keys;
            N=length(keys);
            for i=1:N
                key=keys{i};
                if ~endsWith(key, '.name') && ~strcmp('keys', key)
                    roi=map.get(key);
                    str=MatBasics.XyToString(...
                        RoiUtil.Position(roi));
                    name=map.get([key '.name']);
                    props.setProperty([key '.position'],  str);
                    props.setProperty([key '.name'],  name);
                    props.setProperty([key '.type'],  class(roi));
                end
            end
            fldr=fileparts(file);
            File.mkDir(fldr);
            File.SaveProperties2(file, props);
            if ~asked
                File.OpenFolderWindow(file, '', true);
            end
            ok=true;
        end
        
        function Delete(map, roi)
            [key, name]=LabelLegendRoi.Find(map,roi);
            if ~isempty(key)
                roi=map.remove(key);
                if ~isempty(roi)
                    delete(roi);
                    map.remove([key '.name']);
                end
                fprintf('Deleted ROI, key=%s, name=%s]n', key, name);
            else
                % not identified with key and name yet
                delete(roi);
            end
        end
        
        function Rename(map, roi)
            if ~isvalid(roi)
                return;
            end
            [key,name]=LabelLegendRoi.Find(map,roi);
            if isempty(key)
                isNew=true;
                key=map.get('keys');
                if isempty(key)
                    key=1;
                end
                map.set('keys', key+1);
                key=num2str(0-key);
                map.set(key, roi);
            else
                isNew=false;
            end
            oldClr=RoiUtil.SetColor(roi, 'red');
            [name, cancelled]=inputDlg('Enter ROI name', 'ROI names', name);
            RoiUtil.SetColor(roi, oldClr);
            if ~cancelled
                RoiUtil.SetLabel(roi, name);
                map.set([key '.name'], name);
                if isNew
                    map.set(key, roi);
                end
            end
        end
        
        function [key, name]=Find(map, roi)
            keys=map.map.keys;
            rois=map.map.values;
            N=length(rois);
            for i=1:N
                if isequal(roi, rois{i})
                    key=keys{i};
                    name=map.get([key '.name']);
                    return;
                end
            end
            key='';
            name='';
        end
    end
end