

function Img_RGB = MY_display_color_nothreshold(V2D,bar_value,MyMap)

V2D(isnan(V2D))=0;
map_nothre_reshape=V2D;
defaultMap = MyMap;
nothrebar = [min(bar_value(:)) max(bar_value(:))];
tmap = ones([size(map_nothre_reshape),3]);
map_normalize = (map_nothre_reshape-nothrebar(1))/(nothrebar(2)-nothrebar(1));
map_normalize(map_normalize<=0) = 0.000000001;
map_normalize(map_normalize>=1) = 1;
nan_mask = isnan(map_normalize);
map_normalize(isnan(map_normalize)) = min(bar_value(:));
tmap(:,:,1) = reshape(defaultMap(ceil(map_normalize*64),1)*256,size(map_normalize));
tmap(:,:,2) = reshape(defaultMap(ceil(map_normalize*64),2)*256,size(map_normalize));
tmap(:,:,3) = reshape(defaultMap(ceil(map_normalize*64),3)*256,size(map_normalize));

% mask_index = find(double(nan_mask)==1);
% tmap(mask_index+numel(map_nothre_reshape)*0) = 0;
% tmap(mask_index+numel(map_nothre_reshape)*1) = 0;
% tmap(mask_index+numel(map_nothre_reshape)*2) = 0;

Img_RGB = tmap;
end
