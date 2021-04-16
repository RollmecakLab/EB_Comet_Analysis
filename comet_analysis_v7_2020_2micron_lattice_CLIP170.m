%- this is version 7 of comet analysis in 2020. it is adjusted to only take
%a one channel kymo
%- it was revisited to look at JC's tyr adn detyr data =. Heavily sourced
%from the original comet_analysis.m files. 
% it will

%  read in a direcotry of kymogrpahs, 
% scan anf find the max
% fit a gaussian to the max-1 pixel plus background to center the plots
% use these center positions to align all the comets..
% if it loses the max you will have to click on it in the plot window
% if there are two pixels with equal max brightness it takes the one on the
% left
% output of interest :

% X2 and Ymean are the average comet  
% The cell arrayz All_Aligned_Profiled, All_Ampl. etc.. can be useful to
% plot and individula comet



% experiment setting, must know GR
draw=1; % if draw=1 plot stuff
pixelSize=.077 ; % microns
Growthrate=3; %um/min  -  check these units!!!!!!! needed to transform data later
exposure = .1 ; % camera exposure (s)
secPerFrame = 3; % time between Frames

%- first read in the image you wnat to use, it should be a two color
%kymogrpah where channel TWO represents the +tip tracking protein
%-*** MODIFIED FOR JC data to only need a one channel kymo!!!
folder_name=uigetdir;   %- get dir
 files=dir(fullfile(folder_name,'*.tif')); % - get files
 files= {files.name}; % - get files names
 numFiles=numel(files); % get num of files

for(j=1:numFiles) %- j will loop over each kymo
     close all; % start with no plots
    Tipimage=imread((fullfile(folder_name,'/',files{1,j})),1); % - gets the GFP (channel 2) image
     'tube number' % displays the tube number
     j% shows j
    files{1,j}% displays the name of the current file
    imagedimensions=size(Tipimage);
    h=imagedimensions(1); %- number of lines in the kymoheighty can use to "prune" to a time..
    w=imagedimensions(2); %- width of the kymo for x-axis plotting and shit..
    X=(1:w)'; %- gets the X-values for the plot...

    endtime = 60 ; % the last frame to measure. for now take 180 = 3 mins at 1 fps
    %- now read in the profiles of this particular kymo
    for(k=1:h) % was 1:endtime
        profiles(:,1,k)=X; %- the x-values for the plots
        profiles(:,2,k)=Tipimage(k,1:end)'; %- this reads in all the plots for thi kymo
    end

%- now for each profile, find the peak 

% - start with the user selecting the approximate peak location and go from
% there...





for(m=1:h) % was 1:h
    % which line of the kymo it is
    currentProfile=profiles(:,2,m);
    
    %- NEW PART***********
    % If this is the first line of a new kymo the user must tell it where
    % the comet is. Click using the LEFT mouse button near the comet peak
    
    if(m==1) % for the first line in the kymo
        plot(currentProfile)
        [xpk,ypk]=ginput(1)
        userselectedpk=round(xpk);
        %- now use this to find the true maximum within say 3 pixels either
        %side of where you clicked
        maximumPosition = min(find(currentProfile==max( currentProfile(userselectedpk-3:userselectedpk+3))))
        lastgoodposition=maximumPosition
        
        % - for getting the next max, should be within 5 pixels of last good
        % location - this will filter out errorneous peaks that are bigger but
        % not the comets..may need to adjust the 5 px part...
        
    elseif (lastgoodposition-5 >0 && lastgoodposition+5<w)
        maximumPosition=min(find(currentProfile==max( currentProfile(lastgoodposition : end)))) % only look form last good position to teh end of the kymo on the right
        
        % -here a quick check to make sure things didn't go horrible
        % wrong...
        if(maximumPosition>lastgoodposition+5) % - that is, if the new peak is 10 pixels from the old one we might have a probelm so let the user pick it agian
            % changed to 5 pixels
            
         plot(currentProfile)
        [xpk,ypk]=ginput(1)
        userselectedpk=round(xpk);
        %- now use this to find the true maximum within say 3 pixels either
        %side of where you clicked
        if(userselectedpk<4)
            userselectedpk=4;
        elseif(userselectedpk>w-3)
            userselectedpk=w-6;
        end
        maximumPosition = min(find(currentProfile==max( currentProfile(userselectedpk-3:userselectedpk+3))));        
        end

        lastgoodposition = maximumPosition;
    else
        maximumPosition=min(find(currentProfile==max( currentProfile(lastgoodposition)))); % if cant find use new
        lastgoodposition = 6 % default to thro tomaximumPosition
    end
         
    
    
   % - how far on either side of peak to go . Take 1 pixel on left and ~
   % 1.5 um on right
   numOnRight = 20 ; % - go 20 pixels past the max  ( so ~ 1.5 um in solution) out in to the background
   numOnLeft = 26 ; % 0 num of pixels to go in to the lattice for fun...
   
   if( lastgoodposition-5 >0 && lastgoodposition+5<w)
   maxvalue=max( currentProfile(lastgoodposition-5: lastgoodposition+5));
   else
   maxvalue=max( currentProfile(lastgoodposition)); % if getting close to end just use this as a data dummy holder  
   end
   if(draw==1)
   subplot(2,1,1);
   plot(currentProfile); % - plots the profile and a line at the max
   line([maximumPosition maximumPosition], [0 maxvalue]) ;
   %pause
   end
   % - Get the data to fit now..
   startidx=min(maximumPosition)-numOnLeft;%- in case max has 2 values
   % if there is more than one max take the lower of them
   endidx=min(maximumPosition)+numOnRight; % there should only be 1 in verison 2 with new filter
   if(startidx>0 && endidx<w) % as long as there is enough data on either side then we...
       
       tempProfile(:,1)=profiles(startidx:1:endidx,1,m);
       tempProfile(:,2)=currentProfile(startidx:1:endidx);
       if(draw ==1)
       subplot(2,1,2);
       plot(tempProfile(:,1),tempProfile(:,2));
       end
       
       %    pause
       %- good up until here..not do the fits...
   %----MIGHT HAVE TO ADJUST THE FITTING START PARAMS BASED ON YOUR
   %DATA-----%
   bgEstimate = mean( currentProfile(endidx-15:endidx)); % - an estimate of the background off the comet
  % fitX=currentProfile(maximumPosition-2:maximumPosition+numOnRight,1);
 %  fitY=currentProfile(maximumPosition-2:maximumPosition+numOnRight,2);
 
 % setup the fit ro be from 2 before the max to the end...
   [maxpos,maxidx]=max(tempProfile(:,2));
   if(maxidx>4)
       fitX= tempProfile(maxidx-2:end,1);
       fitY = tempProfile(maxidx-2:end,2);
   else
       maxidx=4;
       fitX= tempProfile(maxidx-2:end,1);
       fitY = tempProfile(maxidx-2:end,2);
   end
%     subplot(3,1,3);
%    plot(fitX,fitY)
      if (numel(fitX)<4)
          g.rsquare=0;
      else
  [f,g]=fit( fitX,fitY,'a1*exp(-((x-b1)/c1)^2)+d1','Startpoint',[max(tempProfile(:,2)), maximumPosition , 1, bgEstimate],'Lower',[ 0 0 0 0]); %-was Startpoint
     end
  % When gaussian fit is good get the params...
  if(g.rsquare>.95)
      FitParams2=coeffvalues(f);
      Amp2=FitParams2(1);
      Center2=FitParams2(2);
      Sigma2=FitParams2(3);
      Background2=FitParams2(4);

      %- for the good ones, let's also get the mean and SD of lattice and
      %background for each 
      
       All_Latticel{m,j}=mean( currentProfile(1: round(Center2)-5 )) % - the average lattice on the left
       All_Lattice_SD{m,j}=std( currentProfile(1: round(Center2)-5 )) % - the average lattice on the left
      
      All_solution{m,j}=mean( currentProfile(round(Center2)+5:end ));
      All_solution_SD{m,j}=std( currentProfile(round(Center2)+5:end ));
      
      
      % for plotting the fit , dont do now..
%       subplot(3,1,3);
%       plot(fitX,fitY)
%       hold on
%       plot(f);   % - plot the temp profile and the fitted data
%       %pause(.2); % can slow down if want to see plotting
%       hold off% so you can look at the plots...
%       %pause(.1);
    
      % - now shift the profiles by the center location....
      
      ShiftedProfile(:,1)=profiles(:,1)-Center2;
      ShiftedProfile(:,2)=profiles(:,2);
      %
      tempShiftProfile(:,1)=tempProfile(:,1)-Center2;
      tempShiftProfile(:,2)=tempProfile(:,2);
      
      %-- good up to here now add the shifted profiled into a cell arrayz..
      All_ALigned_Profiles{m,j}=tempShiftProfile;
      All_Ampl{m,j}=Amp2;
      All_Bg{m,j}=Background2;
      AllSigma2{m,j}=Sigma2;
      clear tempShiftProfile
      clear ShiftedProfile
  end % ends the processing if the profile fits a gausian well..
   %-clear them before next file starts
   end % ends the fitting as long as there is enough data

   end % ends the loop over this file

  clear profiles % to get ready for the next file
end % - ends the loop over all files


%***** NOW WE HAVE ALL THE SHIFTED PROFILEs, WE MUST BIN THEN AVERAGE THEM
%THEN FIT. For the TAIL DROPOFF WE WILL ONLY FIT THE LEFTMOST 20 BINSTHAT
%IS 20PIXELS TO THE LEFT OF THE MAX...

numberOfBins=20; %was numberOneitherSideofPeak;
AllBinnedData=cell(1,numberOfBins+1); %- to hold the binned data now...
NumFileWithTips=(size(All_ALigned_Profiles)); % in case some kymos give no data
for(k=1:NumFileWithTips(2))%-loop over all the files
    numOfProfiles=  numel(All_ALigned_Profiles(:,k)); % - the number of profiles in each column now. bin...
    for(m=1:numOfProfiles) %- for each profile in this file....
        
        if(isempty(All_ALigned_Profiles{m,k})==0) %_only do for non empty ones
            tempprofile(:,1)=All_ALigned_Profiles{m,k}(:,1);
            tempprofile(:,2)=All_ALigned_Profiles{m,k}(:,2); %_ got this holdign the temporary profile now...
            for(j=0:numberOfBins) %_ start with data at zero..
                
                % - find the index of values that go in this bin
                foundidx=find(tempprofile(:,1)>=j & tempprofile(:,1)<j+1);
                %- now put them in place
                for(h=1:numel(foundidx))
                    AllBinnedData{1,j+1}(end+1)=tempprofile(foundidx(h),2);
                end % end loop over all data in this bin
            end % - end loop over all bins
        end
        clear tempprofile %- clear it and start again...
        
    end
    % - end loop over all profiles
end % - end loop over all files


for(q=1:19)
    X2(q)=q.*pixelSize; % - get in to um
    Ymean(q)=mean(AllBinnedData{1,q}); %- FOR A MEAN FILTER
    Ymedian(q)=median(AllBinnedData{1,q}); %-FOR A MEDIAN FILTER BUT WORK ON THIS LATER
end



% - to look at decay LENGTH
hold off;
plot(X2',Ymean')
hold on
scatter(X2',Ymean')
hold on

[f2,g2]=fit( X2',Ymean','a1*exp(-(x/b1))+c1','Startpoint',[Ymean(1)    ,1,1],'Lower',[0 0 0]); %-was Startpoint
Final_FitParams=coeffvalues(f2);
Final_Amp=Final_FitParams(1)
Final_decay=Final_FitParams(2)
Final_bgoffset=Final_FitParams(3)

hold on;
plot(f2)








%_transform from space to time to lookat decay rates...

% 
%  Xtime=X2*pixelSize/Growthrate;
% % - Why do this part??
% 
% %_-now plot and fit
% 
% plot(Xtime',Ymean')
% hold on
% scatter(Xtime',Ymean')
% hold on
% 
% [f2,g2]=fit( Xtime',Ymean','a1*exp(-(x/b1))+c1','Startpoint',[Ymean(1)    ,1,1],'Lower',[0 0 0]); %-was Startpoint
% Final_FitParams=coeffvalues(f2);
% Final_Amp=Final_FitParams(1)
% Final_decay=Final_FitParams(2)
% Final_bgoffset=Final_FitParams(3)
% 
% hold on;
% plot(f2)
% 
%  

 