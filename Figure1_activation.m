clc;clear
close all
WholePath = 'F:\BF_optogentics\';
codepath = 'F:\BF_optogentics\fMRI_code_Trange\';
cd(codepath);
addpath(genpath(codepath));


Type = {'chat';'vglut2';'pv';'som'};
for tloop=1:numel(Type)

Func3D = fullfile(WholePath,'2nd_TwoSampleTTest_results',...
    upper(Type{tloop}),'spmT_0001.nii');

barlim = [-10 -2.2 +2.2 10];
% if tloop>=3;barlim(3)=barlim(3)/1.5;end
Temp3D = fullfile(codepath,'Template_Mouse_v38.nii');
[f1_rgb]=Colormap_3Dviewer_v3(Func3D,Temp3D,barlim,codepath,1,0.70);
[f2_rgb]=Colormap_3Dviewer_v3(Func3D,Temp3D,barlim,codepath,1,0.5);

I = cat(1,f1_rgb.cdata,f2_rgb.cdata);

filename = fullfile(WholePath,'2nd_TwoSampleTTest_results',...
    upper(Type{tloop}),'3D_ACTIVATION.tiff');
imwrite(uint8(I),filename);
close all;


I = f1_rgb.cdata;
filename = fullfile(WholePath,'2nd_TwoSampleTTest_results',...
    upper(Type{tloop}),'3D_spmT_0001.tiff');
imwrite(uint8(I),filename);
close all;

end


Type = {'chat';'vglut';'pv';'sst'};
for tloop=1:numel(Type)

Func3D = fullfile(WholePath,'2nd_TwoSampleTTest_results',...
    upper(Type{tloop}),'spmT_0001.nii');

barlim = [-10 -2.89 +2.89 10];

I0 = spm_read_vols(spm_vol(Func3D));
I = spm_read_vols(spm_vol(fullfile(codepath,'Colormap_3Dviewer','Template_Mouse_X20.nii')));
I0x = imresize3(I0,size(I),'linear');
I0x = flip(flip(I0x,2),3);
I0x = smooth3(I0x,'box',3);

Ip = permute(I,[3 2 1]);
I0xp = permute(I0x,[3 2 1]);

lmask = Ip<10^5;


    
barlim = [-10 10]/1;
if tloop==2;barlim = [-10 10]/1;end
if tloop==1;barlim = [-10 10]/1;end

I0 = spm_read_vols(spm_vol(Func3D));


I = spm_read_vols(spm_vol(fullfile(codepath,'Colormap_3Dviewer','Template_Mouse_X20.nii')));
I0x = imresize3(I0,size(I),'linear');
I0x = flip(flip(I0x,2),3);
I0x = smooth3(I0x,'box',3);
I0x(I0x>-2/1 & I0x<+2/1)=nan;


I0x(I<10000)=nan;
I(:,:,1:10)=[];
I0x(:,:,1:10)=[];
I0x(I0x<barlim(1))=barlim(1);
I0x(I0x>barlim(2))=barlim(2);

I0x(01:10,01:10,01:10)=barlim(1);
I0x(11:20,11:20,11:20)=barlim(2);

for vloop=1:2
    config = struct('CameraPosition',[-4 0 0],...
                    'CameraUpVector',[0,0,1],...
                    'CameraTarget',[0 0 0],...
                    'CameraViewAngle',15,...
                    'BackgroundColor',[0 0 0]+1,...
                    'Renderer','MaximumIntensityProjection',... %%VolumeRendering
                    'Alphamap',(0:255)'/255,...
                    'Lighting',0,...
                    'IsosurfaceColor',[1 1 1]*0.5,...
                    'Isovalue',0.5);
    config.Colormap = [jet(256-20);ones(20,3)];
    if vloop==1
        config.CameraPosition=[0 -3.3 0];
        config.CameraUpVector = [-1 0 0];
    end
    if vloop==2
        config.CameraPosition=[-3.3 0 0];
        config.CameraUpVector = [0 1 0];
    end

    F1=figure;volshow(I,config);
    f1_rgb = getframe(F1);

gdmap = [(0:127)'/127,(0:127)'/127,ones(128,1)];
drmap = [ones(128,1),(0:127)'/127,(0:127)'/127];
defaultMap = [gdmap;flipud(drmap);];

    Map =  jet;%MY_viridis_colormap(256)/256;
    config.Colormap = Map;
    config.Lighting = 0;
    F2=figure;volshow(I0x,config);
    f2_rgb = getframe(F2);

    f1 = f1_rgb.cdata;
    f2 = f2_rgb.cdata;
    f = f1*(1/2)+f2*(1/2);
    if vloop==1;F01=f;end
    if vloop==2;F02=f;end
end
close all
Fa = cat(1,F01,F02);

filename = fullfile(WholePath,'2nd_TwoSampleTTest_results',...
    upper(Type{tloop}),'spmT_0001.tiff');
imwrite(uint8(Fa),filename);
close all;

%}



end