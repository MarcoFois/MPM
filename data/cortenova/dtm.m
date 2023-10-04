clear all
clc

rowstart = 1;
rowend   = 175;

colstart = 1;
colend   = 165;

ndivrows = rowend - rowstart;
ndivcols = colend - colstart;

%[mask_fin,R3] = readgeoraster('mask_out.tif');

immagine = imread ("mask_in_vladi.tif"); %, "PixelRegion", {[rowstart rowend], [colstart, colend]});
%immagine_int = reshape (typecast (immagine(end:-1:1, :), "int16"), size (immagine));

spy(immagine)
find(nonzeros(immagine))

