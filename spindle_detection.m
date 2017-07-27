function varargout=spindle_detection(varargin)
% [timespindles durspindles MINS DENS]=spindle_detection(EDFfile, XMLfile, channels, references, just2, method, fspindle, ecglabel)
% given an EDF and an XML file, detects spindles based on wavelet analysis
% or bandpass analysis
%
% (in the future I will add other features)
%
% MANDATORY INPUTS:
%
% EDFfile
%         string containing file name
% XMLfile
%         string containing file name
% channels
%         cell array containing channels names
% references
%         cell array containing reference channels names
% just2
%         flag (if =1 use only stage 2, otherwise whole sig)
% method
%         flag (if =0 use wavelet method, if =1 bandpass)
%
% OPTIONAL INPUT:
%
% fspindle
%         1x1 double containing central frequency of spindles (default = 13.5)
% ecglabel
%         string containing ecg channel name
%
% OUTPUTS:
%
% timespindles
%              vector of starting point of each spindle, in decimal seconds from
%              beginning of signal
% durspindles
%               vector of duration of each spindle, in decimal seconds
% MINS
%         minutes of clean N2 sleep (vector, 1 value per channel)
% DENS
%         spindle density (vector, 1 value per channel)
%
% Adapted from methods described in Purcell, S.M., 
% Manoach, D.S., Demanuele, C., 
% Cade, B.E., Mariani, S., Cox, R., Panagiotaropoulou, G., 
% Saxena, R., Pan, J.Q., Smoller, J.W. and Redline, S., 
% 2017. Characterizing sleep spindles in 11,630 individuals 
% from the National Sleep Research Resource. 
% Nature Communications, 8.
%
% Developed by Sara Mariani at Brigham and Women's Hospital
% July 2017
% This program is distributed in the hope that it will be useful, but WITHOUT ANY
% WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
% PARTICULAR PURPOSE. 
%
% please report bugs to sara.mariani7@gmail.com

%Set default pararameter values
inputs={'EDFfile','XMLfile','channels','references','just2','method','fspindle','ecglabel'};
outputs={'timespindles','durspindles','MINS','DENS'};

Ninputs=length(inputs);
if nargin>Ninputs
    error('Too many input arguments')
end
if nargin<6
    error('Not enough input arguments')
end

eval([inputs{7} '=13.5;']) %default
eval([inputs{8} '=[];']) %default

for n=1:nargin
    eval([inputs{n} '=varargin{n};'])
end

Noutputs=length(outputs);
if nargout>Noutputs
    error('Too many output arguments')
end

% open signal
if ~isempty(references)
    ch=[channels references];
else
    ch=channels;
end
[header, signalHeader, signalCell] = blockEdfLoad(EDFfile, ch);
% insert here a check for the sampling frequencies;
for kk=1:length(ch)
    SR(kk)=signalHeader(kk).samples_in_record;
end
if numel(unique(SR))~=1
    error('different sampling rates')
end
sttime=header.recording_starttime;
sttime=[str2num(sttime(1:2)) str2num(sttime(4:5)) str2num(sttime(7:8))];
if sttime(1)>12 sttime(1)=sttime(1)-12; end
sttimesec=sttime(1)*3600+sttime(2)*60+sttime(3);

% obtain hypnogram
s=parseXML(XMLfile);
if isfield(s,'Children')
    name=s.Children(5).Name;
    if strcmp(name,'SleepStages')
        hyp=[];
        for k=1:length(s.Children(5).Children)
            hyp(k)=s.Children(5).Children(k).Children.Data;
        end
    else
        return
    end
    if hyp(1)>=48
        hyp=hyp-48;
    end
elseif isfield(s,'C')
    name=s.C(5).Name;
    if strcmp(name,'SleepStages')
        hyp=[];
        for k=1:length(s.C(5).C)
            hyp(k)=s.C(5).C(k).C.Data;
        end
    else
        return
    end
    if hyp(1)>=48
        hyp=hyp-48;
    end
else
    error('Check fields of s')
end
% get events: arousal,movement, artifact
ev=s.C(4).C;
evstart=NaN(size(ev));
evdur=NaN(size(ev));
for jj=1:length(ev)
    c=ev(jj).C;
    name=c(1).C.Data;
    start=c(2).C.Data;
    dur=c(3).C.Data;
    if ~isempty(strfind(name,'artifact')) || ...
            ~isempty(strfind(name,'Arousal')) || ...
            ~isempty(strfind(name,'Movement'))
        evstart(jj)=str2num(start);
        evdur(jj)=str2num(dur);
    end
end
evstart=evstart(~isnan(evstart));
evdur=evdur(~isnan(evdur));
evstart(evstart==0)=1;
eventvec=zeros(length(hyp)*30,1);
for jj=1:length(evstart)
    eventvec(evstart(jj):evstart(jj)+evdur(jj)-1)=1;
end
numepochs=floor(length(eventvec)/30);
eventvec(numepochs*30+1:end)=[];
eventvec=ceil(mean(reshape(eventvec,30,length(eventvec)/30)));

for ch=1:length(channels)
    if ~isempty(references)
        sig=signalCell{ch}-signalCell{ch+length(channels)};
    else
        sig=signalCell{ch};
    end
    fS=signalHeader(ch).samples_in_record;
    t=[1:length(sig)]/fS;
    
    % bandpass filter [0.3-35 Hz]
    d = designfilt('bandpassfir', 'FilterOrder',21,...
        'CutoffFrequency1',0.3,'CutoffFrequency2',35, ...
        'SampleRate', fS);
    sig=filtfilt(d,sig);
    ls=length(sig);
    % remove epochs with arousal, movement or artifact annotation
    epochevent=find(eventvec);
    for jj=1:length(epochevent)
        sig((epochevent(jj)-1)*30*fS+1:epochevent(jj)*30*fS)=NaN;
    end
    sig(ls+1:end)=[];
    % now remove artifacts
    % 30 second epochs
    deltaBand = [1 4];
    betaBand = [15 30];
    numepochs=floor(length(sig)/fS/30);
    delta=zeros(numepochs,1);
    beta=zeros(numepochs,1);
    cleansig=sig;
    for jj=1:numepochs
        win=sig((jj-1)*30*fS+1:jj*30*fS);
        [P,F]=pwelch(win,hamming(4*fS),2*fS,[],fS);
        fbin=F(2)-F(1);
        delta(jj)=sum(P(F>=deltaBand(1) & F<deltaBand(2)))*fbin;
        beta(jj)=sum(P(F>=betaBand(1) & F<betaBand(2)))*fbin;
    end
    deltaavg=smooth(delta,15);
    deltaidx=delta./deltaavg;
    betaavg=smooth(beta,15);
    betaidx=beta./betaavg;
    epochsart=unique([find(deltaidx>2.5);find(betaidx>2)]);
    for jj=1:length(epochsart)
        cleansig((epochsart(jj)-1)*30*fS+1:epochsart(jj)*30*fS)=NaN;
    end
    
    % now remove ECG artifact
    % get ECG
    % can modify later if user enters label of ecg
    
    if ~isempty (ecglabel)
        try
            [~, ecgHeader, ecgsig] = blockEdfLoad(EDFfile,ecglabel);
        catch
            try
                [~, ecgHeader, ecgsig] = blockEdfLoad(EDFfile,{'ECG'});
            catch
                try
                    [~, ecgHeader, ecgsig] = blockEdfLoad(EDFfile,{'EKG'});
                catch
                    [~, ecgHeader, ecgsig] = blockEdfLoad(EDFfile,{'ECG1'});
                end
            end
        end
    else
        try
            [~, ecgHeader, ecgsig] = blockEdfLoad(EDFfile,{'ECG'});
        catch
            try
                [~, ecgHeader, ecgsig] = blockEdfLoad(EDFfile,{'EKG'});
            catch
                [~, ecgHeader, ecgsig] = blockEdfLoad(EDFfile,{'ECG1'});
            end
        end
    end
    cleansigecg=ecgDecont(cleansig,fS,cell2mat(ecgsig),ecgHeader.samples_in_record);
    
    hyps=repmat(hyp,30*fS,1);
    hyps=hyps(:);
    if length(hyps)<length(cleansigecg)
        hyps=[hyps;zeros(length(cleansigecg)-length(hyps),1)];
    end
    
    if just2==1
        % select only EEG during stage 2: resample hyp in EEG samples
        sig2=cleansigecg;
        sig2(hyps~=2)=NaN;
    else
        sig2=cleansigecg;
    end
    
    % remove saturated values
    
    filters=zeros(length(hyp),3);
    satu=zeros(length(hyp),1);
    % statistical filters epoch-by-epoch
    % and remove saturated values
    for jj=1:length(filters)
        win=sig2((jj-1)*30*fS+1:jj*30*fS);
        [a, m, c]=hjorth(win(~isnan(win)));
        filters(jj,:)=[a m c];
        if numel(find(abs(win)==max(abs(sig))))>length(win)/100*5
            satu(jj)=1;
            sig2((jj-1)*30*fS+1:jj*30*fS)=NaN;
        end
    end
    nop=[];
    for jj=1:3 % for each Hjorth parameter
        f=filters(:,jj); %act, mob or comp
        for ii=1:3
            s=nanstd(f);
            nop=[nop;find(f>nanmean(f)+2*s);find(f<nanmean(f)-2*s)];
            f(unique([find(f>nanmean(f)+2*s);find(f<nanmean(f)-2*s)]))=NaN;
        end
    end
    nop=unique(nop);
    for jj=1:length(nop)
        sig2((nop(jj)-1)*30*fS+1:nop(jj)*30*fS)=NaN;
    end
    
    % first plot: signal, stage 2 and artifacts
%     figure
%     ax(1)=subplot(411);
%     plot(t,sig)
%     xlabel('')
%     set(gca,'xtick',[], 'fontsize',15)
%     title(channels(ch))
%     ax(2)=subplot(412);
%     stairs([1:numepochs]*30-30,deltaidx,'r')
%     title('Slow frequency band index')
%     xlabel('')
%     set(gca,'xtick',[], 'fontsize',15)
%     ax(3)=subplot(413);
%     stairs([1:numepochs]*30-30,betaidx,'g')
%     title('Fast frequency band index')
%     xlabel('')
%     set(gca,'xtick',[], 'fontsize',15)
%     ax(4)=subplot(414);
%     plot(t,cleansigecg)
%     hold on
%     plot(t,sig2,'m')
%     legend('Signal after ECG artifact decontamination','signal I use')
%     xlabel('Time (s)')
%     set(gca,'xtick',[], 'fontsize',15)
%     linkaxes(ax,'x')
    if method==1
        % design my wavelet
        fc=13.5;
        n=7;
        ss=n/(2*pi*fc);
        tp=2*ss^2;
        tStart = -4;
        tStop = 4;
        timeVector = linspace(tStart,tStop, (tStop-tStart)*fS );
        psiWavelet = (pi*tp)^(-0.5).*exp(2*1i*pi*fc.*timeVector).*exp(-timeVector.^2/tp);
        
%         figure
%         subplot(121)
%         plot(timeVector,real(psiWavelet), timeVector,imag(psiWavelet));
%         xlabel('Time / Seconds');
%         title(sprintf('Morlet Wavelet, T_{period} = %2.2f, f_c = %2.2f Hz', tp, fc));
        
        input = psiWavelet;
        Nfft = 10 * 2^nextpow2(length(input));
        psd = 20.*log10(fftshift(abs(fft(input,Nfft))));
        freqs = [0:Nfft - 1].*(fS/Nfft);
        freqs(freqs >= fS/2) = freqs(freqs >= fS/2) - fS;
        freqs = fftshift(freqs);
%         subplot(122)
%         plot(freqs, psd);
%         xlabel('Frequency / Hz');
%         title(sprintf('Morlet Wavelet, T_{period} = %2.2f, f_c = %2.2f Hz', tp, fc));
%         grid on
        
        % computing baseline
        t=[1:length(sig2)]/fS;
        tt=t(~isnan(sig2));
        s=sig2(~isnan(sig2));
        hyps2=hyps(~isnan(sig2));
        fsig=conv(s,psiWavelet,'same');
        ssig=smooth(abs(real(fsig)),0.1*fS);
        baseline=mean(ssig(hyps2==2));
        
%         figure
%         ax1(1)=subplot(311);
%         plot([1:length(s)]/fS,s)
%         
%         ax1(2)=subplot(312);
%         plot([1:length(fsig)]/fS,real(fsig))
%         % [cfs,frequencies] = cwt(s,scales,'cmor1-1.5',1/512);
%         % plot([1:length(fsig)]/512,cfs(155,:))
%         % apply moving average of 0.1 s
%         ax1(3)=subplot(313);
%         plot([1:length(ssig)]/fS,real(ssig))
%         linkaxes(ax1,'x')
        
        th=2*baseline;
        cores=ssig>th;
        hold on
        stairs([1:length(ssig)]/fS,cores*500)
        
        % detect consecutive ones between 0.3 and 3 s
        coress=num2str(cores)';
        pattern=num2str(ones(round(fS*0.3),1));
        beg=strfind(coress,pattern'); %start points of each candidate core
        b1=[0;diff(beg')];
        beg(b1==1)=[];
%         plot(beg/fS,ones(size(beg)),'*g')
        
        % store duration of each core
        dur=zeros(size(beg));
        for jj=1:length(beg)
            k=1;
            while(cores(beg(jj)+k)==1 && length(cores)>beg(jj)+k)
                k=k+1;
            end
            dur(jj)=k;
        end
%         plot((beg+dur)/fS,ones(size(beg)),'*m')
        beg(dur>3*fS)=[];
        dur(dur>3*fS)=[];
        
        % now extend the spindles: at least 0.5 s above th2=2*baseline
        th2=1.5*baseline;
        candidates=ssig>th2;
        
        % detect consecutive ones greater than 0.5 s
        candidatess=num2str(candidates)';
        pattern=num2str(ones(round(fS*0.5),1));
        begc=strfind(candidatess,pattern'); %start points of each candidate core
        b1=[0;diff(begc')];
        begc(b1==1)=[];
%         plot(begc/fS,ones(size(begc)),'*c')
        
        % store duration of each extension
        durc=zeros(size(begc));
        for jj=1:length(begc)
            k=1;
            while(candidates(begc(jj)+k)==1 && length(candidates)>begc(jj)+k)
                k=k+1;
            end
            durc(jj)=k;
        end
 %       plot((begc+durc)/fS,ones(size(begc)),'*k')
        
        % now combine the two classifications
        C=zeros(size(cores));
        for jj=1:length(beg)
            C(beg(jj):beg(jj)+dur(jj))=1;
        end
        
        E=zeros(size(cores));
        for jj=1:length(begc)
            E(begc(jj):begc(jj)+durc(jj))=1;
        end
        % intersect
        EC=E&C;
        % adjust duration
        b1=diff(EC');
        begec=find(b1==1);
        
        ind=zeros(size(begec));
        
        for jj=1:length(begec)
            truebeg=find(begc<=begec(jj));
            if ~isempty(truebeg)
            ind(jj)=truebeg(end);
            else ind(jj)=NaN;
            end
        end
        ind(isnan(ind))=[];
        BEG=begc(ind);
        DUR=durc(ind);
        
%         plot(BEG/fS,2*ones(size(BEG)),'>c')
%         plot((BEG+DUR)/fS,2*ones(size(BEG)),'<m')
        
        % now that I found the spindles, next step is merging those too close
        % (closer than 1 sec)
        EPTS=BEG+DUR;
        distance=-EPTS(1:end-1)+BEG(2:end);
        while (any(distance<fS))
            ind=find(distance<fS);
            ind=ind(1);
            if (DUR(ind)+distance(ind)+DUR(ind+1))<3*fS
                BEG(ind+1)=[];
                DUR(ind)=DUR(ind)+distance(ind)+DUR(ind+1);
                DUR(ind+1)=[];
                distance=BEG(2:end)-(BEG(1:end-1)+DUR(1:end-1)); %update distance
            else
                distance(ind)=999; %do not merge
            end
        end
        
%         plot(BEG/fS,2*ones(size(BEG)),'>g')
%         plot((BEG+DUR)/fS,2*ones(size(BEG)),'<r')
%         legend('processed signal','thresholded','core start','core end','estension start', ...
%             'extension end','spindle start','spindle end','spindle start','spindle end')
    else
        t=[1:length(sig2)]/fS;
        tt=t(~isnan(sig2));
        s=sig2(~isnan(sig2));
        hyps2=hyps(~isnan(sig2));
        fs=fS;
        %bandpass filter signal in the 11-15 hz band
        B = fir1(101,[11/fs*2 15/fs*2],'bandpass');
        sigbp=filtfilt(B,1,s);
        
        %calculate the RMS of the bandpass filtered signal
        % with a time resolution of 250 ms
        % using a time window of 250 ms (no overlap)
        numsamples=round(fs/1000*250);
        numwindows=floor(length(sigbp)/numsamples);
        sig=reshape(sigbp,numsamples,numwindows);
        sigrms=rms(sig);
        hyps3=reshape(hyps2,numsamples,numwindows);
        hyps3=round(median(hyps3));
        th=prctile(sigrms(hyps3==2),95);
        cores=sigrms>th;
%         figure
%         ax1(1)=subplot(211);
%         plot([1:length(s)]/fs,s)
%         ax1(2)=subplot(212);
%         plot([1:length(sigrms)]/fs*numsamples,sigrms)
%         hold on
%         stairs([1:length(sigrms)]/fs*numsamples,cores*5)
%         linkaxes(ax1,'x')
%         % detect consecutive ones between 0.3 and 3 s
        beg=strfind(cores,[1 1]); %start points of each candidate core
        b1=[0;diff(beg')];
        beg(b1==1)=[];
        % store duration of each core
        dur=zeros(size(beg));
        for jj=1:length(beg)
            k=1;
            while(cores(beg(jj)+k)==1)
                k=k+1;
            end
            dur(jj)=k;
        end
        beg(dur>3*4)=[];
        dur(dur>3*4)=[];
%         plot(beg/fS*numsamples,ones(size(beg)),'*g')
%         plot((beg+dur)/fS*numsamples,ones(size(beg)),'*m')
        
        % now that I found the spindles, next step is merging those too close
        % (closer than 1 sec)
        BEG=beg;
        DUR=dur;
        EPTS=BEG+DUR;
        distance=-EPTS(1:end-1)+BEG(2:end);
        while (any(distance<fS))
            ind=find(distance<fS);
            ind=ind(1);
            if (DUR(ind)+distance(ind)+DUR(ind+1))<3*fS/32
                BEG(ind+1)=[];
                DUR(ind)=DUR(ind)+distance(ind)+DUR(ind+1);
                DUR(ind+1)=[];
                distance=BEG(2:end)-(BEG(1:end-1)+DUR(1:end-1)); %update distance
            else
                distance(ind)=999; %do not merge
            end
        end
        
%         plot(BEG/fS*32,2*ones(size(BEG)),'>g')
%         plot((BEG+DUR)/fS*32,2*ones(size(BEG)),'<r')
        legend('processed signal','thresholded','spindle start','spindle end')
    end
    beg=tt(BEG);
    DUR=DUR/fS;
    mins=numel(s)/60/fS;
    dens=length(BEG)/mins;
    timespindles{ch}=beg';
    durspindles{ch}=DUR';
    MINS(ch)=mins;
    DENS(ch)=dens;
end

for n=1:nargout
    eval(['varargout{n}=' outputs{n} ';'])
end