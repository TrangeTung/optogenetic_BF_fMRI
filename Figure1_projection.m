clc;clear
close all
WholePath = 'F:\BF_optogentics\';
codepath = 'F:\BF_optogentics\fMRI_code_Trange\';
cd(codepath);
addpath(genpath(codepath));


Type = {'chat';'vglut';'pv';'sst'};
for tloop=1:numel(Type)

Func3D = fullfile(WholePath,'BF_structure_connectivity',...
    ['nR_',Type{tloop},'_with_mask.nii']);

barlim = [ +0 100];

I0 = spm_read_vols(spm_vol(Func3D));
I = spm_read_vols(spm_vol(fullfile(codepath,'Colormap_3Dviewer','Template_Mouse_X20.nii')));
I0x = imresize3(I0,size(I),'linear');
I0x = flip(flip(I0x,2),3);
I0x = smooth3(I0x,'box',5);

Ip = permute(I,[3 2 1]);
I0xp = permute(I0x,[3 2 1]);

lmask = Ip<10^5;

IMG=[];
for sl=57:22:114
    Img_RGB = MY_display_function_map_3D_nothreshold(I0xp,barlim,sl);
    Img_RGB = flip(Img_RGB,1);
    Lx = flip(lmask(:,:,sl)',2);
    SE = strel('disk',2); Lx = imdilate(Lx,SE);
    Lx = 1-Lx;
    Lx3 = repmat(Lx,[1 1 3]);
    Img_RGB = Img_RGB.*double(Lx3);
    Img_RGB(Lx3~=1)=255;
    IMG = cat(1,IMG,Img_RGB);
end
IMG(:,end-4:end,:)=255;
figure;imshow(uint8(IMG),[])

filename = fullfile(WholePath,'BF_structure_connectivity',...
    ['nR_',Type{tloop},'_activations.tiff']);
imwrite(uint8(IMG),filename);
close all;



end