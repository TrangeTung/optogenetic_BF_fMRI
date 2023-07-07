clc;clear
codepath = 'F:\BF_optogentics\fMRI_code_Trange\';
cd(codepath);
addpath(genpath(codepath));
WholePath = 'F:\BF_optogentics\\';

Parc = fullfile(WholePath,'voxelwise_connectome','01_connectome',...
    'mapping_voronoi_to_100um','final_tesselation.nii.gz');
ihdr = spm_vol(Parc);
Pres = spm_read_vols(ihdr);
Pv = fmask(Pres,ones(size(Pres)));


load(fullfile(WholePath,'voxelwise_connectome','01_connectome',...
    'full_connectome_percolation_thr.mat'));
V = -log10(full_connectome);
V(full_connectome==0)=0;
V = max(V(:))-V;
V(full_connectome==0)=0;
%%%%%%%%%%%%%%%%%%%%
V0 = V(1:size(V,1)/2,1:size(V,2)/2)+V(1:size(V,1)/2,1+size(V,2)/2:end);
V=V0; clear V0 full_connectome
%%%%%%%%%%%%%%%%%%%%


dest = fullfile(WholePath,'DimReduc_SC');
Type = {'chat';'pv';'sst';'vglut';};
DMx = [];
for tloop=1:numel(Type)
    
    filename = fullfile(WholePath,'BF_structure_connectivity',...
        ['R_',Type{tloop},'_with_mask.nii']);
    Func_Img_3D = spm_read_vols(spm_vol(filename));
    Trace = ( Func_Img_3D + flip(Func_Img_3D,1) )/2;
    Trace = flip(Trace,3);
    Vs = zeros(size(V,1),1);
    for loop=1:size(V,1)
        lmask = Pres==loop;
        Vs(loop) = nanmean(fmask(Trace,lmask))/255*100;
    end
    DMx = [DMx,Vs];
end
DMx(isnan(DMx))=0;
DMx(DMx~=0)=1;


y = V; 
TypeName={'Nature';'Montaged';'Chat';'PV';'SOM';'Vglut2'};

for tl = 2%:6
    
    if tl==2;DM=[y.*repmat(DMx(:,1),1,size(V,1));y.*repmat(DMx(:,2),1,size(V,1));...
            y.*repmat(DMx(:,3),1,size(V,1));y.*repmat(DMx(:,4),1,size(V,1))];end
    

    if tl==1;DM=y;end
    if tl==3;DM=y.*repmat(DMx(:,1),1,size(V,1));end
    if tl==4;DM=y.*repmat(DMx(:,2),1,size(V,1));end
    if tl==5;DM=y.*repmat(DMx(:,3),1,size(V,1));end
    if tl==6;DM=y.*repmat(DMx(:,4),1,size(V,1));end
    
    DM0=DM';
    R2 = nan(64,1);
    Wall = cell(64,1);
    Hall = cell(64,1);
    
    %
    for st=1:64
        
        fprintf(['BF No. ',num2str(st),'\n']);
        
        DM0=DM';
        DM0(:,sum(DM0,1)==0)=[];
        [W,H,errs,loss] = nmf_euc(DM0', st);
        
        ySC = W*H;

        Vec = DM0(:)'; VecP = ySC(:)';
        mu = mean(VecP,2).*mean(Vec,2) - mean(Vec.*VecP,2);
        md = mean(VecP,2).^2 - mean(VecP.^2,2);
        b = mean(Vec,2)-mu./md.*mean(VecP,2);
        Vecf = mu./md.*VecP + b;
        SSR = sum((Vecf-mean(Vec,2)).^2,2);
        SST = sum((Vec-mean(Vec,2)).^2,2);
        R2(st,1) = SSR./SST;
        Wall{st,1}=W;
        Hall{st,1}=H;
        
    end
    mkdir(fullfile(WholePath,'DimReduc_SC_log_voxwise_new'));
    cd(fullfile(WholePath,'DimReduc_SC_log_voxwise_new'));
    
    save(['Celltype_R2_6W_',TypeName{tl},'.mat'],'R2');
    save(['Celltype_W_6W_',TypeName{tl},'.mat'],'Wall','-v7.3');
    save(['Celltype_H_6W_',TypeName{tl},'.mat'],'Hall','-v7.3');
    %}
    
    cd(fullfile(WholePath,'DimReduc_SC_log_voxwise_new'));
    
    load(['Celltype_R2_6W_',TypeName{tl},'.mat']);
    load(['Celltype_W_6W_',TypeName{tl},'.mat']);
    load(['Celltype_H_6W_',TypeName{tl},'.mat']);
    
    Waligned = cell(6,1);
    Haligned = cell(6,1);
    R2aligned = cell(6,1);
    for w=1:64
        
        if w~=1
            ref = Wall{w-1};
            mov = Wall{w};
            C = abs(corr(ref,mov));
            ALLpin = zeros(w-1,1);
            for lp=1:w-1
                pin = find(C(lp,:)==max(C(lp,:)));
                ALLpin(lp)=pin;
                C(:,pin)=0;
            end
            ALLpin(end+1)=setdiff(1:w,ALLpin);
            Wall{w} = mov(:,ALLpin);
            Hall{w} = Hall{w}(ALLpin,:);
        end
        
        Wx = Wall(w,:); Hx = Hall(w,:);
        R2collec = zeros(w,1);
        loop=1;
        Wxx = Wx{loop};
        Hxx = Hx{loop};

        R2 = zeros([w,1]);
        for xl=1:w
            ySC = Wxx(:,xl)*Hxx(xl,:);
            
            DM0=DM';
            DM0(:,sum(DM0,1)==0)=[];

            Vec = DM0(:)'; VecP = ySC(:)';
            mu = mean(VecP,2).*mean(Vec,2) - mean(Vec.*VecP,2);
            md = mean(VecP,2).^2 - mean(VecP.^2,2);
            b = mean(Vec,2)-mu./md.*mean(VecP,2);
            Vecf = mu./md.*VecP + b;
            SSR = sum((Vecf-mean(Vec,2)).^2,2);
            SST = sum((Vec-mean(Vec,2)).^2,2);
            R2(xl) = SSR./SST;
        end
        [sR2,s] = sort(R2,'descend');
        R2collec(:,loop) = sR2;

        Wxy = Wxx(:,s);
        Wx{loop} = Wxy;
        Hxy = Hxx(s,:);
        Hx{loop} = Hxy;
        
        Wr = cell2mat(Wx);
        Wmean = zeros(size(Wxx));
        for xl=1:w
            Wmean(:,xl) = mean(Wr(:,xl:w:end),2);
        end
        Waligned{w} = Wmean;
        Hr = cat(1,Hx{:});
        Hmean = zeros(size(Hxx));
        for xl=1:w
            Hmean(xl,:) = mean(Hr(xl:w:end,:),1);
        end
        Haligned{w} = Hmean;
        R2aligned{w} = mean(R2collec,2);
    end
    
    save(['Celltype_Waligned_',TypeName{tl},'.mat'],'Waligned');
    save(['Celltype_Haligned_',TypeName{tl},'.mat'],'Haligned');
    save(['Celltype_R2aligned_',TypeName{tl},'.mat'],'R2aligned');
    
    
end


%
TypeName={'Nature';'Montaged';'Chat';'PV';'SOM';'Vglut2'};
dest = fullfile(WholePath,'DimReduc_SC_log_voxwise');
ihdr = spm_vol(fullfile(dest,'Template_Mouse_X20_I2.nii'));
for tl = 1:6
    cd(dest);
    load(['Celltype_Waligned_',TypeName{tl},'.mat']);
    
    
    for wl=1:10    
        
        filepath = fullfile(dest,TypeName{tl},['N',num2str(wl)]);
        mkdir(filepath);
        
    	W = Haligned{wl}';
        W=[W;W];
        st=wl;
        W4d = zeros([size(Pres),st]);
        Wv = fmask(W4d,ones(size(Pres)));
        for lp=1:size(W,1)
            for ix=1:st
            Wv(Pv==lp,ix) = W(lp,ix);
            end
        end
        W4d = funmask(Wv,ones(size(Pres)));
        for ix=1:st
           hdr(ix)=ihdr;
           hdr(ix).fname = fullfile(filepath,'NMF.nii');
           hdr(ix).n=[ix,1];
           hdr(ix).dt=[16,0];
        end
        spm_write_vol_4D(hdr,W4d);
        


        barlim = [-4 -1 1 4];
        Temp3D = fullfile(codepath,'Template_Mouse_v38.nii');

        ihdr = spm_vol(fullfile(filepath,'NMF.nii'));
        I = spm_read_vols(ihdr);
        I0=[];
        for lp=1:size(I,4)
            I0(:,:,:,lp) = imresize3(I(:,:,:,lp),[114 80 38]);
            I0(:,:,:,lp)=flip(I0(:,:,:,lp),3);
            hdr = spm_vol(Temp3D);
            ihdr(lp).mat=hdr.mat;
            ihdr(lp).dim=hdr.dim;
            ihdr(lp).fname =fullfile(filepath,['cut_NMF.nii']);
        end
        spm_write_vol_4D(ihdr,I0);




        for ix=1:wl
            template = Parc;
            spmT_file = fullfile(filepath,['NMF.nii,',num2str(ix)]);
            Colormap(   'statfile',spmT_file,...
                'bgmfile',template,...
                'slice',102:-7:9,...
                'bar_value',[-5 -2 2 5]*1,...
                'dest',filepath,...
                'mapname',['NMF_C',num2str(ix)],...
                'cluster',10);
            
            
            
            Func3D = fullfile(filepath,['cut_NMF.nii']);
            %loc={[0.45,0.3];[1,0.75];[1,0.5];[0.75,0.5];[0.5,0.3];[1,0.75]};
            loc = [1;0.8;0.6;0.4];
            [f1_rgb]=Colormap_3Dviewer_v3(Func3D,Temp3D,barlim,codepath,ix,loc(1));
            [f2_rgb]=Colormap_3Dviewer_v3(Func3D,Temp3D,barlim,codepath,ix,loc(2));
            [f3_rgb]=Colormap_3Dviewer_v3(Func3D,Temp3D,barlim,codepath,ix,loc(3));
            [f4_rgb]=Colormap_3Dviewer_v3(Func3D,Temp3D,barlim,codepath,ix,loc(4));
            
            
            I1 = cat(1,f2_rgb.cdata);
            I2 = cat(1,f3_rgb.cdata,f4_rgb.cdata);
            I0 = cat(1,I1,I2);
            
            filename = fullfile(filepath,['3D_2_ICA_',num2str(ix),'.tiff']);
            imwrite(uint8(I0),filename);
            close all
        end
        
        
        
    
    end
end
%}
