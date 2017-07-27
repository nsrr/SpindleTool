function allfilespindles(datafolder,resultfolder, ...
    analysissignals, referencesignals, referencemethod, ...
    only2, ecgname, fspin, start, xmlsuffix, spinmethod)
%
% computes spindle detection in all the selected
% files, saves results in spreadsheets
% Developed by Sara Mariani at Brigham and Women's Hospital
% July 2017
% This program is distributed in the hope that it will be useful, but WITHOUT ANY
% WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
% PARTICULAR PURPOSE. 
%
% please report bugs to sara.mariani7@gmail.com
if referencemethod==1
    referencesignals=repmat(referencesignals,length(analysissignals),1);
end
toplabels=cell(1,length(analysissignals)*2);
toplabels(1:2:end)=analysissignals;
labels=repmat([{'Start_time (s)'},{'Duration (s)'}],1,length(analysissignals));
files=dir([datafolder './*.edf']);
files=files(start:end);
M=zeros(length(files),length(analysissignals));
D=zeros(length(files),length(analysissignals));
for jj=1:length(files)
    EDFfile=[datafolder '\' files(jj).name];
    names{jj}=files(jj).name(1:end-4);
    XMLfile=[EDFfile(1:end-3) xmlsuffix];
    [timespindles durspindles MINS DENS]=spindle_detection(EDFfile, ...
        XMLfile, analysissignals, referencesignals, only2, spinmethod, ...
        fspin, ecgname);
    % create spreadsheet for this individual subject
     totalPowerFn = ([resultfolder files(jj).name(1:end-4) '_spindletimes.xlsx']);
     % each lead will have a different numer of spindles
     % take the highest number
     le=zeros(1,length(analysissignals));
     for ii=1:length(analysissignals)
         le(ii)=length(timespindles{ii});
     end
     totalSummary=cell(max(le),length(analysissignals)*2);
     for ii=1:length(analysissignals)
     totalSummary(1:le(ii),(ii-1)*2+1)=cellstr(num2str(timespindles{ii}));
     totalSummary(1:le(ii),(ii-1)*2+2)=cellstr(num2str(durspindles{ii}));
     end
     totalSummary=[toplabels;labels;totalSummary];
     xlswrite(totalPowerFn, totalSummary);
     M(jj,:)=MINS;
     D(jj,:)=DENS;
     save MD M D
end
% create summary spreadsheet for all subjects
totalPowerFn =([resultfolder '\spindlesummary.xlsx']);
lab=cell(1,length(analysissignals)*2+1);
lab(1)={'Subject'};
lab(2:2:end)=analysissignals;
lab2=cell(1,length(analysissignals)*2+1);
lab2(2:2:end)={'Minutes of usable signal'};
lab2(3:2:end)={'Spindle density'};
for ii=1:length(analysissignals)
MM(:,(ii-1)*2+1)=cellstr(num2str(M(:,ii)));
MM(:,(ii-1)*2+2)=cellstr(num2str(D(:,ii)));
end
mymat=[names MM];
mymat=[lab;lab2;mymat];
xlswrite(totalPowerFn, mymat);