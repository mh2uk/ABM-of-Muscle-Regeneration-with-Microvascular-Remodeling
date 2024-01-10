%% Process histology image for ABM

clear all; close all; clc;

%% Read in images
img = imread('exampleImage.png'); % input histology image here
img = img(1:3:end,1:3:end,:); % down sample image
% read in binary image
bw = imread('exampleImageBinary.png'); % input binarized histology image here
bw = bw(1:3:end,1:3:end,:);
bw = imbinarize(bw(:,:,1));

%% Identify ECM, fiber, and membranes 
fiber_membranes = bwperim(bw);
image_w_membranes = double(bw) + double(fiber_membranes);
imshow(image_w_membranes./2)

% currently in image_w_membrane
% ECM = 0; fibers = 1; membranes = 2
[row, col] = find(image_w_membranes == 0);
ECM = [row, col, ones(length(row), 1)];

[row, col] = find(image_w_membranes == 1);
Fibers = [row, col, ones(length(row), 1).*2];

[row, col] = find(image_w_membranes == 2);
Fiber_membranes = [row, col, ones(length(row), 1).*3];


all_codes = [ECM; Fibers; Fiber_membranes];
all_codes = [all_codes(:,2), all_codes(:,1), all_codes(:,3)];
% all_codes(:,2) = (all_codes(:,2) - 120).*-1;
all_codes = sortrows(all_codes);

pointsize = 5;
scatter(all_codes(:,1), all_codes(:,2), pointsize, all_codes(:,3), 'filled');
axis equal

%% Save csv
csvwrite('exampleCSV.csv', all_codes)
