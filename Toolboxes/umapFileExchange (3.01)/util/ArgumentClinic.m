classdef ArgumentClinic < handle
    properties(SetAccess=private)
        fig;
        app;
        H;
        argsObjs=cell(1,3);
        jDlgClinics=cell(1,3);
    end
    methods
        function this=ArgumentClinic
            this.app=BasicMap.Global;
            this.fig=Gui.NewFigure(true, 'off');
            fp=get(this.fig, 'OuterPosition');
            if this.app.highDef
                fp(3)=fp(3)*1.015;
                fp(4)=fp(4)*.680;
            else
                fp(3)=fp(3)*1.15;
                fp(4)=fp(4)*.80;
            end
            set(this.fig, 'OuterPosition', fp);
            movegui(this.fig, 'north');
            this.fig.Name='Stanford University''s Herzenberg pipelines...';
            Gui.JavaComponentAt(ArgumentClinic.TopPanel(this.app,...
                @(idx)selectClinic(this, idx)), 'north', this.fig);
            this.fig.Visible='on'; 
        end
        
        function selectClinic(this, idx)
            if ~isempty(this.jDlgClinics{idx})
                this.jDlgClinics{idx}.setVisible(true);
                this.jDlgClinics{idx}.requestFocus;
            elseif ~isempty(this.argsObjs{idx})
                argsObj=this.argsObjs{idx};
            else
                if idx==2
                    argsObj=SuhMatch.GetArgsWithMetaInfo();
                    where='south++';
                    forRunningWhat='for running QFMatch/QF-tree';
                    ttl='subset characterization';
                elseif idx==3
                    [file, example]=SuhEpp.EliverArgs;
                    argsObj=SuhEpp.GetArgsWithMetaInfo(file, example{:});
                    where='south east++';
                    forRunningWhat='for running EPP';
                    ttl='unsupervised subset identification';
                else
                    example={'label_column', 'end', ...
                        'match_scenarios', 2};
                    argsObj=UmapUtil.GetArgsWithMetaInfo(...
                        'eliverLabeled.csv',  example{:});
                    forRunningWhat='for running UMAP/UST';
                    ttl='supervised subset identification';
                    where='south west++';
                end
                
                jd=javax.swing.JDialog;
                jd.getContentPane.add(Gui.FlowPanelCenter(15,11,...
                    argsObj.getArgumentClinicPanel(forRunningWhat)));
                jd.setTitle(['SUH ' ttl]);
                jd.pack;
                this.jDlgClinics{idx}=jd;
                SuhWindow.Follow(this.jDlgClinics{idx}, ...
                    this.fig, where, true);
                jd.setVisible(true)
            end
        end
    end
    
    methods(Static)
        function pnl=TopPanel(app, fnc)
            if nargin<2
                app=BasicMap.Global;
                if nargin<1
                    fnc=[];
                end
            end
            [~,leonard]=Gui.ImageLabel(Html.Wrap(Html.ImgXy(...
                'Leonard.png', [], .5, false, false, app)), [], ...
                'See about our lab''s founder', @seeLen);
            
            arthur=Gui.Label(['<html>'...
                Html.ImgXy('facs.gif',[], .94, false,false, app) ...
                '</html>']);
            pipeline=Gui.Label(['<html>'...
                Html.ImgXy('pipeline.png',[], .94, false,false, app) ...
                '</html>']);
            wayne=html('wayneMoore2.png', .39, ...
                '<u>unsupervised</u> identification', 'exhaustive projection pursuit (EPP)', 'Wayne Moore');
            connor=html('connor.png',.15, '<u>supervised</u> identification', 'parameter reduction (UMAP/UST)', 'Connor Meehan');
            orlova= html('darya.png',.2, ...
                'characterization', 'QFMatch/QF-tree', 'Darya Orlova');
            items={ connor, orlova, wayne, orlova};
            jRadio=Radio.PanelAndCallback(@resolve, true, items);
            normalFont=javax.swing.UIManager.getFont('Label.font');
            font=java.awt.Font(normalFont.getName, ...
                normalFont.ITALIC, normalFont.getSize+3);
            pnlLen=Gui.BorderPanel([], 0, 8, ...
                'North', leonard, 'South', ...
                Html.WrapSmallBold('Len Herzenberg'));
            pnlSouthWest=Gui.BorderPanel;
            pnlSouthWest.add(pnlLen, 'West');
            pnlSouthWest.add(Gui.Panel(arthur), 'East');
            if app.highDef
                pnl=Gui.Panel( Gui.BorderPanel([], 2, 15, ...
                    'North', pipeline, 'South', ...
                    pnlSouthWest), Gui.BorderPanel([], 0, 2, 'North',...
                    Html.WrapHr(['<font color="blue"><b>'...
                    'Choose a data subsetting pipeline</b></font>']),...
                    'Center', jRadio));
            else
                pnl=Gui.Panel( Gui.BorderPanel([], 2, 15, ...
                    'North', pipeline, 'South', ...
                    pnlSouthWest), jRadio);
                Gui.SetTitledBorder('Choose a data subsetting pipeline', pnl, font);
            end
            function str=html(img, scale, words, via, provider)
                img=Html.ImgXy(img,[], scale, ...
                    false,false,app);
                str=Html.Wrap(['<table cellspacing="5"><tr><td>' ...
                    img '</td><td><font color="blue"><b>Subset '...
                    words '</b></font><br>via <i>', via ...
                    '</i><br>' Html.WrapBoldSmall([' ' provider])...
                    '</td></tr></table>']);
            end
            
            function seeLen(h,e)
                web('https://www.pnas.org/content/110/52/20848', '-browser');
            end
            
            function resolve(h,e)
                item=char(e.getActionCommand);
                idx=StringArray.IndexOf(items, item);
                if isempty(fnc)
                    fprintf('Index %d chosen\n', idx);
                else
                    feval(fnc, idx);
                end
            end
        end
    end
end