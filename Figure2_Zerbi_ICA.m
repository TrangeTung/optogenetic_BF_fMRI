clc;clear
codepath = 'F:\BF_optogentics\fMRI_code_Trange\';
cd(codepath);
addpath(genpath(codepath));
WholePath = 'F:\BF_optogentics\\';
NII_v213 = fullfile(codepath,'Label_Mouse_213_v38.nii');
xhdr = spm_vol(NII_v213);
Labels_v38 = spm_read_vols(xhdr);


filename = fullfile(WholePath,'Zerbi','ABI_DR_zerbi15_200um_reorient.nii.gz');
ihdr = spm_vol(filename);
I = spm_read_vols(ihdr);

Ia = zeros([114 80 38 size(I,4)]);
clear rhdr
for il=1:size(I,4)
    Ir = imresize3(I(:,:,:,il),[114 80 38]);
    Ir = flip(flip(Ir,3),2);
    Ia(:,:,:,il)=Ir;
    rhdr(il) = xhdr;
    rhdr(il).fname = fullfile(WholePath,'Zerbi','ABI_DR_zerbi15.nii');
    rhdr(il).dt=[16,0];
    rhdr(il).n = [il;1];
end
spm_write_vol_4D(rhdr,Ia);



ICAname = {'Somatosensory';'Sensory(VC+AUD+Mot)';'DMN';'HPF';'BasalGanglia';'Olfactory'};
ICAs = {[03:06];[02 07 08];[09:10];[11:12];[13:16];[01 17 18]};
Ib = zeros([114 80 38 numel(ICAs)]);
clear bhdr;
for icl=1:numel(ICAs)
    Ib(:,:,:,icl)=mean(Ia(:,:,:,ICAs{icl}),4);
    bhdr(icl) = xhdr;
    bhdr(icl).fname = fullfile(WholePath,'Zerbi','ICA_zerbi15.nii');
    bhdr(icl).dt=[16,0];
    bhdr(icl).n = [icl;1];
end
spm_write_vol_4D(bhdr,Ib);



Func3D = fullfile(WholePath,'Zerbi','ICA_zerbi15.nii');
barlim = [-10 -3 +3 10];
Temp3D = fullfile(codepath,'Template_Mouse_v38.nii');
for lp=3%1:numel(ICAs)
    [f1_rgb]=Colormap_3Dviewer_v3(Func3D,Temp3D,barlim,codepath,lp,1.0);
    [f2_rgb]=Colormap_3Dviewer_v3(Func3D,Temp3D,barlim,codepath,lp,0.6);
    
    if lp==4
        [f1_rgb]=Colormap_3Dviewer_v3(Func3D,Temp3D,barlim,codepath,lp,0.4);
        [f2_rgb]=Colormap_3Dviewer_v3(Func3D,Temp3D,barlim,codepath,lp,0.3);
    end
    if lp==3
        [f1_rgb]=Colormap_3Dviewer_v3(Func3D,Temp3D,barlim,codepath,lp,0.8);
        [f2_rgb]=Colormap_3Dviewer_v3(Func3D,Temp3D,barlim,codepath,lp,0.6);
    end
    if lp==5
        [f1_rgb]=Colormap_3Dviewer_v3(Func3D,Temp3D,barlim,codepath,lp,0.7);
        [f2_rgb]=Colormap_3Dviewer_v3(Func3D,Temp3D,barlim,codepath,lp,0.6);
    end

    I = cat(2,f1_rgb.cdata,f2_rgb.cdata);
    filename = fullfile(WholePath,'Zerbi',['3D_2_ICA_',num2str(lp),'.tiff']);
    imwrite(uint8(I),filename);
    close all
    
    
    template = [codepath,'\meanT2.nii'];
    spmT_file = [Func3D,',',num2str(lp)];
    lmask = [codepath,'\lmask_Mouse_v38.nii'];
    Colormap(   'statfile',spmT_file,...
        'bgmfile',template,...
        'slice',9:2:35,...
        'bar_value',barlim,...
        'dest',fullfile(WholePath,'Zerbi'),...
        'mapname',['xxx2D_ICA_',num2str(lp)],...
        'denoi_profile',lmask,...
        'cluster',10);%,...
    
end







lmask = Labels_v38>1;
Iv = fmask(Ib,lmask);
Im = zeros(size(Iv,1),1);

for ix=1:size(Iv,1)
    I0_ = Iv(ix,:);
    I0_(I0_<5)=0;
    I0 = find(I0_==max(I0_));
    if numel(I0)==1;Im(ix,1)=I0;end
end
I3 = funmask(Im,lmask);
hdr = xhdr; 
hdr.fname = fullfile(WholePath,'Zerbi','ICA_5.nii');
hdr.dt=[16,0];
spm_write_vol(hdr,I3);


