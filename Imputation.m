clear ; close all
scrsz = get(0,'ScreenSize') ;



%folder = 'AgeWise_HTW/T2_Post-CBTI_N=30/' ;
%ggg = dir([folder 'AP1-T2-*.xlsx']) ;

folder = 'SIR_HTW/' ;
ggg = dir([folder 'SIR-T1-*.xlsx']) ;

PLOT = 0 ;

%% step 0: get data
if PLOT
    figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)]) ;
end

fprintf('Step 0: get data\n')

for jj = 1: length(ggg)
    disp(jj)
    [aa, bb, ~] = xlsread([folder ggg(jj).name]) ;
    temp0 = aa(:, end) ;
    HH = hour(aa(:,1)) ;
    MM = minute(aa(:,1)) ;
    if strcmp(folder, 'SIR_HTW/') && jj == 20 
        HH = hour(bb(4:4+length(temp0)-1,1)) ;
        MM = minute(bb(4:4+length(temp0)-1,1)) ;
    end

    if jj == 30 
        HH(1:3) = [] ; MM(1:3) = [] ; temp0(1:3) = [] ;
        aa(1:3,:) = [] ;
    end

    JUMP = 0 ; 
    for kk = 2: length(HH)
        if HH(kk) == 0 && MM(kk) == 0 
            JUMP = JUMP + 24 ;
        end

        if strcmp(folder, 'AgeWise_HTW/T2_Post-CBTI_N=30/') && jj == 19 && kk == 2000 
            JUMP = JUMP + 24 ;
        end

        HH(kk) = HH(kk) + JUMP ;
    end

    % time0 = the min-index of sampled data coming from the xls file
    % Note that there are some jumps (might come from device changes)
    time0 = (HH-HH(1))*60 + MM ;

    % time = the real time with the unit min
    % this is uniform time samples WITHOUT jumps
    time = [1:(time0(end)-time0(1)+1)] + time0(1) - 1 ;

    % temp1 = inital data with potential jumps in time/samples in the xls file
    temp1 = nan(size(time)) ;
    temp1(time0-time0(1)+1) = temp0 ;
    idx0 = find(isnan(temp1)) ;
    idx = find(~isnan(temp1)) ;
    temp1 = interp1(idx, temp1(idx), 1:length(temp1), 'linear') ;


    % Rule1: < 35.5 deg or >39
    idx1 = find(temp1 <= 35.5 | temp1 >= 39) ;
    if strcmp(folder, 'AgeWise_HTW/T2_Post-CBTI_N=30/') && jj == 23
        idx1 = find(temp1 <= 34) ;
    end

    % Rule2: >0.4 deg per 3 mins    
    [dtemp, ~] = derivative(time, temp1) ;
    idx2n = find(dtemp < -0.4/3) ;
    idx2p = find(dtemp > 0.4/3) ;
    idx2 = union(idx2n, idx2p) ;


   

    % extend each found segment by L=10, R=20
    idx = union(idx1, idx2) ;
    %if ~isempty(idx3) ; idx = union(idx, idx3) ; end
    tmpQ = ones(size(temp1)) ; tmpQ(idx) = nan ;
    tmpi = isnan(tmpQ) ;
    [M, start] = regexp(sprintf('%i', [tmpi']), '1+', 'match'); 
    M = cellfun(@length, M);
    for kk = 1: length(M)
        if M(kk) >= 1
            idx = [max(1, start(kk)-10):start(kk)-1 idx ...
                start(kk)+M(kk):min(length(temp1), start(kk)+M(kk)+20)] ;
        end
    end

    idx = union(idx, idx0) ;


    temp = temp1 ;
    temp(idx) = nan ;



    idx1 = isnan(temp) ;
    [M, start] = regexp(sprintf('%i', [idx1']), '1+', 'match'); 
    M = cellfun(@length, M);
    if start(1) == 1
        Starttime = [hour(aa(M(1)+1,1)) minute(aa(M(1)+1,1)) second(aa(M(1)+1,1))] ;
        temp(1:M(1)) = [] ;
        time(1:M(1)) = [] ;
    else
        Starttime = [hour(aa(1,1)) minute(aa(1,1)) second(aa(1,1))] ;
    end

   
    

    %{
    subplot(3,1,mod(jj-1,3)+1)
    hold off ; 
    if start(1)==1
        plot(temp1(M(1)+1:end), 'color', [.7 .7 .7]) ;
    else
        plot(temp1, 'color', [.7 .7 .7]) ; 
    end
    hold on ; plot(temp, 'r', 'linewidth', 2) ;
    if jj~=23 ; axis([-inf inf 35 39]) ; else; axis([-inf inf 34 38]) ; end 
    set(gca,'fontsize', 16) ; 
    filename = ggg(jj).name ;
    filename(filename == '_') = ' ' ;
    title(['Case ' num2str(jj) ': ' filename])
    plot(35.5*ones(size(temp1)), 'b') ; 
    if mod(jj,3)==0 || jj == length(ggg)
        export_fig(['RemovedSegments' num2str(jj)], '-transparent', '-pdf') ; 
        %pause
    end
    %}

    save(ggg(jj).name(1:end-5), 'temp', 'time', 'Starttime', 'temp0', 'time0')
end






%% step 1: get #nights and start recording time

fprintf('Step 1: get #nights and start recording time\n')

if strcmp(folder, 'SIR_HTW/')
    [aa,bb,cc] = xlsread('SIR_vis_scored_sleep_forHTW.xlsx','Visually Scored Sleep_75') ;
end

if strcmp(folder, 'AgeWise_HTW/T2_Post-CBTI_N=30/')
    [aa,bb,cc] = xlsread('AgeWise_GNT_GMT_forHTW.xlsx','PROJ1_11b_EX2_PSG') ;
end

NoNight = zeros(length(ggg), 1) ;
for jj = 1: length(ggg)
    if strcmp(folder, 'AgeWise_HTW/T2_Post-CBTI_N=30/')
        tmp = find(ggg(jj).name(8:end)=='_') ;
    elseif strcmp(folder, 'SIR_HTW/')
        tmp = find(ggg(jj).name(8:end)=='-') ;
    end
    ID = ggg(jj).name(8:8-1+tmp(1)-1) ;
    NoNight(jj) = length(find(aa(:,1)==str2num(ID))) ;
end



% this is obtained from the raw data (the absolute time is saved there)
% NOT from the xls file.
RecordStart = -1*ones(length(aa(:,1)), 4) ;
for jj = 1: length(aa(:,1))
    disp(jj)

    CASEidx = [] ;
    for kk = 1: length(ggg)
        if strcmp(folder, 'AgeWise_HTW/T2_Post-CBTI_N=30/')
            tmp = find(ggg(kk).name(8:end)=='_') ;
        elseif strcmp(folder, 'SIR_HTW/')
            tmp = find(ggg(kk).name(8:end)=='-') ;
        end
        ID = ggg(kk).name(8:8-1+tmp(1)-1) ;
        if str2double(ID) == aa(jj,1)
            CASEidx = [CASEidx kk] ;
        end
    end

    if length(CASEidx) > 1 || isempty(CASEidx)
        fprintf('ERROR\n')
        disp(CASEidx)
    else        
        load([ggg(CASEidx).name(1:end-5) '.mat']) ;
        RecordStart(jj, :) = [Starttime(1:3) CASEidx] ;
    end
end

if strcmp(folder, 'SIR_HTW/')
    save('RecordStartSIR', 'RecordStart', 'NoNight')
end

if strcmp(folder, 'AgeWise_HTW/T2_Post-CBTI_N=30/')
    save('RecordStartPT1', 'RecordStart', 'NoNight')
end





%% step 1-2: construct the library
% this library is for imputation


fprintf('Step 1-2: get library\n')


LIB = {} ;
LIBmean = zeros(length(ggg), 1) ;
N = 1 ;
for jj = 1: length(ggg)
    load([ggg(jj).name(1:end-5) '.mat']) ;
    x0 = temp ;
    idx1 = ~isnan(x0) ;
    
    [M, start] = regexp(sprintf('%i', [idx1']), '1+', 'match'); 
    M = cellfun(@length, M);
    for kk = 1: length(M)
        LIB{N} = x0(start(kk): start(kk)+M(kk)-1) - LIBmean(jj) ;
        LIBmean(N) = mean(x0(idx1)) ;
        N = N + 1 ;
    end
end

load LIB2
for jj = 1: length(LIB2)
    LIB2{jj} = LIB2{jj}' ;
end
LIB = [LIB LIB2] ;
LIBmean = [LIBmean; LIBmean2] ;

if strcmp(folder, 'SIR_HTW/')
    save('LIBSIR', 'LIB', 'LIBmean')
end

if strcmp(folder, 'AgeWise_HTW/T2_Post-CBTI_N=30/')
    save('LIBPT1', 'LIB', 'LIBmean')
end




%% step 2: impute the siganl

fprintf(['Step 2: impute the signal ' num2str(length(ggg)) 'cases\n'])


if strcmp(folder, 'SIR_HTW/')
    load('LIBSIR')
end

if strcmp(folder, 'AgeWise_HTW/T2_Post-CBTI_N=30/')
    load('LIBPT1')
end

if PLOT
    figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)]) ;
end

DATA = {} ;
FILENAME = {} ;

for jj = 1: length(ggg)
    disp(jj)
    load([ggg(jj).name(1:end-5) '.mat']) ;
    filename = ggg(jj).name ;
    filename(filename == '_') = ' ' ;
    x0 = temp ;
    
    idx1 = find(~isnan(x0)) ;
    idx2 = find(isnan(x0)) ;

    % sample period = min
    % use 1 hour for bdry 
    [Imp] = impute(x0, 60, LIB, LIBmean) ;

    if PLOT
        subplot(3,1,mod(jj-1,3)+1)
        hold off ;
        plot(time/60, Imp, 'r', 'linewidth', 2) ; hold on ;
        plot(time/60, x0, 'color', [.7 .7 .7], 'linewidth', 2) ;
        set(gca, 'fontsize', 18) ; ylabel('Temp (C)') ;
        axis tight ; xlabel('Time (h)')
        title(['Case ' num2str(jj) ': ' filename])
        if mod(jj,3)==0 || jj == length(ggg)
            export_fig(['Imputed' num2str(jj)], '-transparent', '-pdf') ;
            pause(1)
        end
        %export_fig(['Case' num2str(jj)], '-transparent', '-pdf') ;
    end

    DATA{jj} = [Imp] ;
    FILENAME{jj} = filename ;



    if strcmp(folder, 'SIR_HTW/')
        csvwrite(['SIR_' filename 'IMP.csv'], Imp) ;
        csvwrite(['SIR_' filename 'Starttime.csv'], Starttime) ;
    end

    if strcmp(folder, 'AgeWise_HTW/T2_Post-CBTI_N=30/')
        csvwrite(['AgeWise_' filename 'IMP.csv'], Imp) ;
        csvwrite(['AgeWise_' filename 'Starttime.csv'], Starttime) ;
    end

end

if strcmp(folder, 'SIR_HTW/')
    save('ImputedDataSIR', 'DATA', 'FILENAME') ;
end

if strcmp(folder, 'AgeWise_HTW/T2_Post-CBTI_N=30/')
    save('ImputedDataPT1', 'DATA', 'FILENAME') ;
end
