function DecayCorrelation2
% This code written 28 December 2017 for the purpose of determing the
% exponential best fit values for the correlation data 
 tic

%Version 1.0 basic functionality, bug fixes
%Version 2.0 custom guess fitting implemented
%Version 3.0 custom guess algorithm perfected, multistart key provides
%global solution with greater consistency 

%Current as of 26 January 2017
 

close all; clearvars; clc
%extract the wells for which data has been recorded
load('EGF(E6)Binned_Correlation_Data.mat','welldatastore');
u = ~cellfun('isempty',welldatastore); %#ok<USENS>
[wells,~,~] = find(u(:,1)==1);
data = cell(size(welldatastore,1),1);

storeALL = nan(numel(wells)*size(welldatastore,2),6);
matstore = cell(size(welldatastore,1),size(welldatastore,2));

mastercounter = 1;

for well = wells'
    tempdatastore = nan(size(welldatastore,2),2);

    for frame = 1:size(welldatastore,2)
        tempwelldata = welldatastore{well,frame}; %#ok<IDISVAR>
%         figure;

        %insert check for NaN values
        tempwelldata(~any(~isnan(tempwelldata), 2),:)=[];
        
        
%         F = @(x,xdata)x(1)*exp(x(2)*xdata);
%         x0 = [0.7,-0.02];
%         [x,resnorm,~,exitflag,output] = lsqcurvefit(F,x0,tempwelldata(:,1),tempwelldata(:,2));
%         tempdatastore(frame,:) = x;
%         param = x;
%         meanres = sum((mean(tempwelldata(:,2))-tempwelldata(:,2)).^2);
%         rsquare = 1-(resnorm/meanres);

        F = @(x,xdata)x(1)*exp(x(2)*xdata);
%         lb = [0.00001 -5];
%         ub = [50 -0.0000001];
        x0 = [0.5 -0.01];
        
        problem = createOptimProblem('lsqcurvefit','x0',x0,'objective',F,...
            'xdata',tempwelldata(:,1),'ydata',tempwelldata(:,2));
        
        ms = MultiStart('PlotFcns',@gsplotbestf);
        [xmulti,errormulti] = run(ms,problem,2000);
        meanres = sum((mean(tempwelldata(:,2))-tempwelldata(:,2)).^2);
        rsquare = 1-(errormulti/meanres);
        x = xmulti;
        
%         [x,resnorm,~,exitflag,output] = lsqcurvefit(F,x0,tempwelldata(:,1),tempwelldata(:,2));
        tempdatastore(frame,:) = x;
        param = x;

       
          tempstore = [well frame param(1) param(2) -1*log(2)/param(2) rsquare];
          matstore{well,frame} = tempstore;
          storeALL(mastercounter,:) = tempstore;
          mastercounter = mastercounter + 1;
            

    end
    

    data{well,1} = tempdatastore;
    disp(strcat('Well',num2str(well),'is now complete'))
end
toc
save('EGF(E6)FitDataTrial2000.mat','storeALL','matstore')
end