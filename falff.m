clear all
close all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%STM-fMRI

AllVolume=readnifti('masklow2mni.nii');
AllVolume=AllVolume(:,:,1:25,:);
MaskData=readnifti('masklow2mni_brain_mask.nii');
MaskData=MaskData(:,:,1:25);
[nDim1 nDim2 nDim3 nDimTimePoints]=size(AllVolume);
BrainSize = [nDim1 nDim2 nDim3];
AllVolume=reshape(AllVolume,[],nDimTimePoints)';
MaskDataOneDim=reshape(MaskData,1,[]);
AllVolume=AllVolume(:,find(MaskDataOneDim));
ASamplePeriod=3;
 nDimTimePoints = size(AllVolume,1);
 sampleFreq = 1/ASamplePeriod;
sampleLength = nDimTimePoints;
paddedLength = 2^nextpow2(sampleLength);

CUTNUMBER = 10;

LowCutoff=0.01;
HighCutoff=0.08;
if (LowCutoff >= sampleFreq/2) % All high included
    idx_LowCutoff = paddedLength/2 + 1;
else % high cut off, such as freq > 0.01 Hz
    idx_LowCutoff = ceil(LowCutoff * paddedLength * ASamplePeriod + 1);
    % Change from round to ceil: idx_LowCutoff = round(LowCutoff *paddedLength *ASamplePeriod + 1);
end
if (HighCutoff>=sampleFreq/2)||(HighCutoff==0) % All low pass
    idx_HighCutoff = paddedLength/2 + 1;
else % Low pass, such as freq < 0.08 Hz
    idx_HighCutoff = fix(HighCutoff *paddedLength *ASamplePeriod + 1);
    % Change from round to fix: idx_HighCutoff	=round(HighCutoff *paddedLength *ASamplePeriod + 1);
end

% Detrend before fft as did in the previous alff.m
%AllVolume=detrend(AllVolume);
% Cut to be friendly with the RAM Memory
SegmentLength = ceil(size(AllVolume,2) / CUTNUMBER);
for iCut=1:CUTNUMBER
    if iCut~=CUTNUMBER
        Segment = (iCut-1)*SegmentLength+1 : iCut*SegmentLength;
    else
        Segment = (iCut-1)*SegmentLength+1 : size(AllVolume,2);
    end
    AllVolume(:,Segment) = detrend(AllVolume(:,Segment));
end
% Zero Padding
AllVolume = [AllVolume;zeros(paddedLength -sampleLength,size(AllVolume,2))]; %padded with zero

%AllVolume = 2*abs(fft(AllVolume))/sampleLength;
% Cut to be friendly with the RAM Memory
SegmentLength = ceil(size(AllVolume,2) / CUTNUMBER);
for iCut=1:CUTNUMBER
    if iCut~=CUTNUMBER
        Segment = (iCut-1)*SegmentLength+1 : iCut*SegmentLength;
    else
        Segment = (iCut-1)*SegmentLength+1 : size(AllVolume,2);
    end
    AllVolume(:,Segment) = 2*abs(fft(AllVolume(:,Segment)))/sampleLength;
end

ALFF_2D = mean(AllVolume(idx_LowCutoff:idx_HighCutoff,:));

% Get the 3D brain back
ALFFBrain = zeros(size(MaskDataOneDim));
ALFFBrain(1,find(MaskDataOneDim)) = ALFF_2D;
ALFFBrain = reshape(ALFFBrain,nDim1, nDim2, nDim3);


% Also generate fALFF
% fALFF_2D = sum(AllVolume(idx_LowCutoff:idx_HighCutoff,:)) ./ sum(AllVolume(2:(paddedLength/2 + 1),:));
fALFF_2D = sum(AllVolume(idx_LowCutoff:idx_HighCutoff,:),1) ./ sum(AllVolume(2:(paddedLength/2 + 1),:),1); %YAN Chao-Gan, 171218. In case there is only one point
fALFF_2D(~isfinite(fALFF_2D))=0;

% Get the 3D brain back
fALFFBrain = zeros(size(MaskDataOneDim));
fALFFBrain(1,find(MaskDataOneDim)) = fALFF_2D;
fALFFBrain = reshape(fALFFBrain,nDim1, nDim2, nDim3);

save fALFFBrain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%NSTM-fMRI

clear all
AllVolume=readnifti('masknstm2std.nii');
AllVolume=AllVolume(:,:,1:25,:);
MaskData=readnifti('masknstm2std_brain_mask.nii');
MaskData=MaskData(:,:,1:25);
[nDim1 nDim2 nDim3 nDimTimePoints]=size(AllVolume);
BrainSize = [nDim1 nDim2 nDim3];
AllVolume=reshape(AllVolume,[],nDimTimePoints)';
MaskDataOneDim=reshape(MaskData,1,[]);
AllVolume=AllVolume(:,find(MaskDataOneDim));
ASamplePeriod=3;
 nDimTimePoints = size(AllVolume,1);
 sampleFreq = 1/ASamplePeriod;
sampleLength = nDimTimePoints;
paddedLength = 2^nextpow2(sampleLength);

CUTNUMBER = 10;

LowCutoff=0.01;
HighCutoff=0.08;
if (LowCutoff >= sampleFreq/2) % All high included
    idx_LowCutoff = paddedLength/2 + 1;
else % high cut off, such as freq > 0.01 Hz
    idx_LowCutoff = ceil(LowCutoff * paddedLength * ASamplePeriod + 1);
    % Change from round to ceil: idx_LowCutoff = round(LowCutoff *paddedLength *ASamplePeriod + 1);
end
if (HighCutoff>=sampleFreq/2)||(HighCutoff==0) % All low pass
    idx_HighCutoff = paddedLength/2 + 1;
else % Low pass, such as freq < 0.08 Hz
    idx_HighCutoff = fix(HighCutoff *paddedLength *ASamplePeriod + 1);
    % Change from round to fix: idx_HighCutoff	=round(HighCutoff *paddedLength *ASamplePeriod + 1);
end

% Detrend before fft as did in the previous alff.m
%AllVolume=detrend(AllVolume);
% Cut to be friendly with the RAM Memory
SegmentLength = ceil(size(AllVolume,2) / CUTNUMBER);
for iCut=1:CUTNUMBER
    if iCut~=CUTNUMBER
        Segment = (iCut-1)*SegmentLength+1 : iCut*SegmentLength;
    else
        Segment = (iCut-1)*SegmentLength+1 : size(AllVolume,2);
    end
    AllVolume(:,Segment) = detrend(AllVolume(:,Segment));
end
% Zero Padding
AllVolume = [AllVolume;zeros(paddedLength -sampleLength,size(AllVolume,2))]; %padded with zero

%AllVolume = 2*abs(fft(AllVolume))/sampleLength;
% Cut to be friendly with the RAM Memory
SegmentLength = ceil(size(AllVolume,2) / CUTNUMBER);
for iCut=1:CUTNUMBER
    if iCut~=CUTNUMBER
        Segment = (iCut-1)*SegmentLength+1 : iCut*SegmentLength;
    else
        Segment = (iCut-1)*SegmentLength+1 : size(AllVolume,2);
    end
    AllVolume(:,Segment) = 2*abs(fft(AllVolume(:,Segment)))/sampleLength;
end

ALFF_2D = mean(AllVolume(idx_LowCutoff:idx_HighCutoff,:));

% Get the 3D brain back
ALFFBrainnstm = zeros(size(MaskDataOneDim));
ALFFBrainnstm(1,find(MaskDataOneDim)) = ALFF_2D;
ALFFBrainnstm = reshape(ALFFBrainnstm,nDim1, nDim2, nDim3);


% Also generate fALFF
% fALFF_2D = sum(AllVolume(idx_LowCutoff:idx_HighCutoff,:)) ./ sum(AllVolume(2:(paddedLength/2 + 1),:));
fALFF_2D = sum(AllVolume(idx_LowCutoff:idx_HighCutoff,:),1) ./ sum(AllVolume(2:(paddedLength/2 + 1),:),1); %YAN Chao-Gan, 171218. In case there is only one point
fALFF_2D(~isfinite(fALFF_2D))=0;

% Get the 3D brain back
fALFFBrainnstm = zeros(size(MaskDataOneDim));
fALFFBrainnstm(1,find(MaskDataOneDim)) = fALFF_2D;
fALFFBrainnstm = reshape(fALFFBrainnstm,nDim1, nDim2, nDim3);

save fALFFBrainnstm


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% show the results
close all
brain=readnifti('masklow2mni_brain.nii');
load('fALFFBrainnstm');
load('fALFFBrain')
total=fALFFBrain-fALFFBrainnstm;
% total1=total>=0.53;
j=0;
for i=1:21
    j=j+1;
figure(j)
imshow(mat2gray(brain(:,:,i))+total1(:,:,i));

end

