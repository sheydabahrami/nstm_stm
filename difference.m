clc
clear all
close all
%obtaining data
datastm=readnifti('masklow2mni.nii');
brain=readnifti('masklow2mni_brain.nii');
% mask=readnifti('maskbinDLPFC');
[x,y,z,t]=size(datastm);
fid=fopen('maskbinDLPFC.img');
mask=fread(fid,[x,y*z]);
mask=reshape(mask,x,y,z);
% % ImP are the voxels in brain
ImP=find(mask==128);
data=[];
for i=1:t
    data1(:,:,:)=datastm(:,:,:,i);
    data11=data1(ImP);
    data=[data data11];
end

[num dim]=size(data);

meandata=mean(data);
figure
plot(meandata)
grid on

DLPFC=meandata;
save DLPFC

corrstm=zeros(size(datastm,1),size(datastm,2),size(datastm,3));
for i=1:size(datastm,1)
    for j=1:size(datastm,2)
        for k=1:size(datastm,3)
            temp=zeros(1,size(datastm,4));
            for p=1:size(datastm,4)
                temp(p)=datastm(i,j,k,p);
            end
            corrstm(i,j,k)=correlation(temp,DLPFC);
          
        end
    end
end
corrstm=corrstm(:,:,1:25);
save corrstm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
%obtaining data
datanstm=readnifti('masknstm2std.nii');
[x,y,z,t]=size(datanstm);
% mask=readnifti('maskbinDLPFC');
fid=fopen('maskbinDLPFC.img');
mask=fread(fid,[x,y*z]);
mask=reshape(mask,x,y,z);
mask=mask(:,:,1:81);

% ImP are the voxels in brain
ImP=find(mask==128);
data=[];
for i=1:t
    data1(:,:,:)=datanstm(:,:,:,i);
    data11=data1(ImP);
    data=[data data11];
end

[num dim]=size(data);

meandata=mean(data);
figure
plot(meandata)
grid on

DLPFC=meandata;
corrnstm=zeros(size(datanstm,1),size(datanstm,2),size(datanstm,3));
for i=1:size(datanstm,1)
    for j=1:size(datanstm,2)
        for k=1:size(datanstm,3)
            temp=zeros(1,size(datanstm,4));
            for p=1:size(datanstm,4)
                temp(p)=datanstm(i,j,k,p);
            end
            corrnstm(i,j,k)=correlation(temp,DLPFC);
          
        end
    end
end
corrnstm=corrnstm(:,:,1:25);
save corrnstm
load ('corrnstm')
NSTM=corrnstm;

load ('corrstm')
STM=corrstm;


j=0;
dif=STM-NSTM;
for i=1:25
    j=j+1;
figure(j)
imshow(dif(:,:,i));

end
j=0;
difnew=dif>=0.5;
for i=1:25
    j=j+1;
figure(j)
imshow(mat2gray(brain(:,:,i))+difnew(:,:,i));

end

