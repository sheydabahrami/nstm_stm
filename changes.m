clc
clear all
close all
load('CORRSTMNEW');
stm=reshape(CORRSTMNEW,91*109*25,1);
load('CORRNSTMNEW');
nstm=reshape(CORRNSTMNEW,91*109*25,1);
figure
plot(nstm,stm,'.')
xlabel('correlation of brain voxels without stimulus')
ylabel('correlation of brain voxels with stimulus')
 grid on
