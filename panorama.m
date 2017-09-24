
clear all;
imagefiles = dir('imgs/2_*.jpg');
nfiles = length(imagefiles)    % Numero de imagens
tforms(nfiles) = projective2d(eye(3));

for i=1:nfiles
    I{i} = rgb2gray(imread(strcat('imgs/', imagefiles(i).name)));
    %imshow(I{i})
end
for i=int32(nfiles/2):2    % Centro para à esquerda
    tforms(i-1) = imgDiff(I{i}, I{i-1})
    tforms(i-1).T = tforms(i-1).T * tforms(i).T
end
for i=int32(nfiles/2):(nfiles-1)   % Centro para à direita
    tforms(i) = imgDiff(I{i}, I{i+1})
    tforms(i+1).T = tforms(i+1).T * tforms(i).T
end

Tinv = invert(tforms(int32(nfiles/2)));
imageSize = size(I{1});
for i = 1:numel(tforms)
    tforms(i).T = tforms(i).T * Tinv.T;
end

for i = 1:numel(tforms)
    [xlim(i,:), ylim(i,:)] = outputLimits(tforms(i), [1 imageSize(2)], [1 imageSize(1)]);
end

% Find the minimum and maximum output limits
xMin = min([1; xlim(:)]);
xMax = max([imageSize(2); xlim(:)]);

yMin = min([1; ylim(:)]);
yMax = max([imageSize(1); ylim(:)]);

%Width and height of panorama.
width  = round(xMax - xMin);
height = round(yMax - yMin);

% Initialize the "empty" panorama.
panoram = []

blender = vision.AlphaBlender('Operation', 'Binary mask', ...
    'MaskSource', 'Input port');

% Create a 2-D spatial reference object defining the size of the panorama.
xLimits = [xMin xMax];
yLimits = [yMin yMax];
panoramaView = imref2d([height width], xLimits, yLimits);

% Create the panorama.
for i = 1:nfiles

    Im = I{i};

    % Transform Im into the panorama.
    warpedImage = imwarp(Im, tforms(i), 'OutputView', panoramaView);

    % Generate a binary mask.
    mask = imwarp(true(size(I,1),size(I,2)), tforms(i), 'OutputView', panoramaView);

    % Overlay the warpedImage onto the panorama.
    panoram = step(blender, panoram, warpedImage, mask);
end

figure
imshow(panoram)

%%%%%%%%%%%%%%%
% Função que retorna matriz de transformação
%%%%%%%%%%%%%%%
function [tform] = imgDiff(orig, dist)
    F1 = detectSURFFeatures(orig);
    F2 = detectSURFFeatures(dist);
    %imshow(I); hold on;
    %plot(F.selectStrongest(10));
    [featuresOriginal,validPtsOriginal] = extractFeatures(orig, F1);
    [featuresDistorted,validPtsDistorted] = extractFeatures(dist, F2);
    index_pairs = matchFeatures(featuresOriginal, featuresDistorted);
    matchedPtsOriginal  = validPtsOriginal(index_pairs(:,1));
    matchedPtsDistorted = validPtsDistorted(index_pairs(:,2));
    %figure; showMatchedFeatures(orig, dist, matchedPtsOriginal, matchedPtsDistorted);
    [tform,inlierPtsDistorted,inlierPtsOriginal] = estimateGeometricTransform(matchedPtsDistorted,matchedPtsOriginal,...
    'projective', 'Confidence', 90.0, 'MaxNumTrials', 3000);
end
