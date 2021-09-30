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
AllVolume=AllVolume(:,find(MaskDataOneDim))';
load DLPFC
corrstm=zeros(74139,1);
for i=1:size(AllVolume,1)
            temp=zeros(1,110);
            for j=1:110
                temp(j)=AllVolume(i,j);
            end
            corrstm(i)=correlation(temp,DLPFC);
          
end
CORRSTMNEW = zeros(size(MaskDataOneDim));
CORRSTMNEW(1,find(MaskDataOneDim)) = corrstm;
CORRSTMNEW = reshape(CORRSTMNEW,nDim1, nDim2, nDim3);
save CORRSTMNEW

clear all
AllVolume=readnifti('masknstm2std.nii');
AllVolume=AllVolume(:,:,1:25,:);
MaskData=readnifti('masknstm2std_brain_mask.nii');
MaskData=MaskData(:,:,1:25);
[nDim1 nDim2 nDim3 nDimTimePoints]=size(AllVolume);
BrainSize = [nDim1 nDim2 nDim3];
AllVolume=reshape(AllVolume,[],nDimTimePoints)';
MaskDataOneDim=reshape(MaskData,1,[]);
AllVolume=AllVolume(:,find(MaskDataOneDim))';
load DLPFC
corrnstm=zeros(66188,1);
for i=1:size(AllVolume,1)
            temp=zeros(1,110);
            for j=1:110
                temp(j)=AllVolume(i,j);
            end
            corrnstm(i)=correlation(temp,DLPFC);
          
end
CORRNSTMNEW = zeros(size(MaskDataOneDim));
CORRNSTMNEW(1,find(MaskDataOneDim)) = corrnstm;
CORRNSTMNEW = reshape(CORRNSTMNEW,nDim1, nDim2, nDim3);
save CORRNSTMNEW
load CORRSTMNEW
load CORRNSTMNEW
DIFF=CORRSTMNEW-CORRNSTMNEW;
j=0;
for i=1:25
    j=j+1;
figure(j)
imshow(DIFF(:,:,i));

end

close all
j=0;
brain=readnifti('masklow2mni_brain.nii');
difnew=DIFF>=0.68;
for i=1:25
    j=j+1;
figure(j)
imshow(mat2gray(brain(:,:,i))+difnew(:,:,i));

end
