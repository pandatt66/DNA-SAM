function sepdist = sep_dist_timeblcks(sam_eqtime, sepeff=0.9973, nsamples=5e5, blcklen=12.0, nsubblcks=1, base = {'dCMP', 'dTMP', 'dGMP', 'dAMP'})

    % Description
    % Calculate distances to separate velocity & flight time distributions of different dNMP pairs. 
    % The maximum separation distance is the required channel length. 
    % First calculates the distribution of mean velocities over some block length using moving block bootstrap, 
    % then determines the number (N) of mean velocities required to separate the overall mean velocity distributions to some separation efficiency. 
    % This is an approximation to running simulations N times longer than the current ones. 
    % Each set of new data with the same length as the existing data would have a different mean velocity whose distribution is estimated using the existing data.
    % A secondary result of the calculation is an estimate of the standard error for the mean velocity.
    %
    % Inputs
    % sam_eqtime: Equilibration time.
    % sepeff: Separation efficiency. The final distributions will have an overlap area of approximately sepeff.
    % nsamples: Number of bootstrap samples for estimating distribution of mean velocity. 
    %           This number needs to be large to avoid variability in the results if calculated multiple times.
    % blcklen: Length of blocks in ns. 
    %          Should be chosen to be longer than autocorrelation time for the velocity estimated where the result of g_analyze starts to level off 
    %          or when the velocity autocorrelation function is nearly zero. 
    %          Alternatively, when the block length is large enough then the resulting channel length will stop changing as a function of block length. 
    %          If the block length is too short, then the channel length will increase with increasing block length.
    % nsubblcks: Number of sub blocks to divide the main block length of blcklen into. Old option. Should be kept at 1.
    % base: Names of dNMP directories. Sub directory names are specified in the read_dNMP_xpos function and will have to be changed for different cases.
    %
    % Outputs
    % Writes mean velocity and standard error to file.
    % Writes separation distances to file.
    % Writes channel length to file.
    % Plots estimated final time of flight distributions.
    % Saves estimated final time of flight distributions data to file.

    %number of dNMPs
    nbases = length(base);
    timecol=1;
    for ibase = 1:nbases
		filenm = [base{ibase} '_1_coord_fixed.dat'];
		t_x_dNMP{ibase} = dlmread(filenm)(:, 1:2);
		ind = t_x_dNMP{ibase}(:,timecol) >= sam_eqtime;
    t_x_dNMP{ibase} = t_x_dNMP{ibase}(ind,:);
		t_x_dNMP{ibase}(:, 2) = -t_x_dNMP{ibase}(:, 2);
	end

    %adjust block lengths to use all data
    bl = blcklen;
    blcklen = zeros(nbases,1);
    nblcks = zeros(nbases,1);
    data_duration = zeros(nbases,1);
    for ibase = 1:nbases
        data_duration(ibase) = t_x_dNMP{ibase}(end,1) - t_x_dNMP{ibase}(1,1);
        nblcks(ibase) = floor((data_duration(ibase)-1e-10)/bl);
        blcklen(ibase) = data_duration(ibase)/nblcks(ibase);
    end
        
    %divide blocks into sub-blocks
    blcklen_short = blcklen/nsubblcks;
    
    for ibase = 1:nbases

        %length of data
        len_xpos = size(t_x_dNMP{ibase})(1);

        %maximum number of data points for one block
        dt_min = min(diff(t_x_dNMP{ibase}(:,1)));
        max_points = ceil(blcklen(ibase)/dt_min);

        %replicate data to deal with cases where time of flight goes past the end of the original data (wrap around)
        xpos_tmp = t_x_dNMP{ibase};
        time_shift = t_x_dNMP{ibase}(end,1) - xpos_tmp(1,1);
        xpos_tmp(:,1) = xpos_tmp(:,1) + time_shift;
        pos_shift = t_x_dNMP{ibase}(end,2) - xpos_tmp(1,2);
        xpos_tmp(:,2) = xpos_tmp(:,2) + pos_shift;
        xpos_tmp = xpos_tmp(2:end,:);
        t_x_dNMP{ibase} = [t_x_dNMP{ibase}; xpos_tmp(1:max_points,:)];
        
        %generate all possible blocks and calculate their velocities
        pos2 = interp1(t_x_dNMP{ibase}(:,1),t_x_dNMP{ibase}(:,2),t_x_dNMP{ibase}(1:len_xpos,1) + blcklen_short(ibase));
        vel_allblocks = (pos2 - t_x_dNMP{ibase}(1:len_xpos,2))/blcklen_short(ibase);
        clear pos2
        
        %choose random block velocities
        nblcks_total = nblcks(ibase)*nsamples; %number of blocks total = number of blocks in the data * number of bootstrap samples
        rndind = randi(len_xpos,nblcks_total,1);
        vel{ibase} = vel_allblocks(rndind);
        
    end

    %reshape velocity data into matrices
    for ibase = 1:nbases
        vel{ibase} = reshape(vel{ibase},length(vel{ibase})/nsamples,nsamples);
    end
    
    % calculate mean of mean velocities, standard deviation of mean velocities, duration of simulations, distance traveled by dNMPs in simulations
    m = zeros(1,nbases);
    s = zeros(1,nbases);
    dt = zeros(1,nbases);
    dx = zeros(1,nbases);
    %~v = zeros(1e7,nbases);
    
    for ibase = 1:nbases
        m_dist{ibase} = mean(vel{ibase})';
        m(ibase) = mean(m_dist{ibase})
        s(ibase) = std(m_dist{ibase});
        dt(ibase) = t_x_dNMP{ibase}(end,1) - t_x_dNMP{ibase}(1,1) - blcklen_short(ibase);
        dx(ibase) = m(ibase)*dt(ibase);
        %~v(ibase,:) = normrnd(m(ibase),s(ibase),1,1e7);
    end

    outdata = [];
    for ibase = 1:nbases
        outdata = [outdata m_dist{ibase}];
    end
    save('mean_vel_dist.dat', 'outdata', '-ascii')
        
    % save mean velocity and standard error
    outfile = ['vel_mean_se_' num2str(sepeff) '.dat'];
    outdata = [m s];
    save(outfile,'outdata','-ascii')
    
    % calculate separation distances for dNMP pairs
    z = abs(norminv((1-sepeff)/2));
    zsq = z^2;

    ns = zeros(nbases,nbases);
    sepdist = zeros(nbases,nbases);

    for ibase = 1:nbases-1
        for jbase = ibase+1:nbases

            ns(ibase,jbase) = zsq*(s(ibase) + s(jbase)*sqrt((m(ibase)*dt(ibase))/(m(jbase)*dt(jbase))))^2/(m(ibase)-m(jbase))^2;
            ns(jbase,ibase) = ns(ibase,jbase)*(m(ibase)*dt(ibase))/(m(jbase)*dt(jbase));
            sepdist(ibase,jbase) = ns(ibase,jbase)*m(ibase)*dt(ibase);
            sepdist(jbase,ibase) = sepdist(ibase,jbase);

        end
    end
    
    % channel length
    channel_length = max(max(sepdist));
    
    % save separation distances & channel length in microns to file
    outfile = ['sepdist_' num2str(sepeff) '.dat'];
    outdata = sepdist/1000;
    save(outfile,'outdata','-ascii');
    
    outfile = ['channel_length_' num2str(sepeff) '.dat'];
    outdata = channel_length/1000;
    save(outfile,'channel_length','-ascii');
    
    % standard deviation for final distributions0
    s_final = zeros(nbases,1);
    for ibase = 1:nbases
        s_final(ibase) = s(ibase)/sqrt(channel_length/dx(ibase));
    end
    
    % analysis time per dNMP
    dist_edges = zeros(nbases,2);
    for ibase = 1:nbases
        dist_edges(ibase,:) = [channel_length/(m(ibase) + z*s_final(ibase)), channel_length/(m(ibase) - z*s_final(ibase))];
    end    
    t_analysis = max(dist_edges(:,2)) - min(dist_edges(:,1));
    
    % save analysis time per dNMP
    outfile = ['tanalysis_' num2str(sepeff) '.dat'];
    outdata = [t_analysis min(dist_edges(:,1))];
    save(outfile,'outdata','-ascii');
    
    % time of flight probability distribution functions
    for ibase = 1:nbases
        x{ibase} = linspace(channel_length/(m(ibase) + 4*s_final(ibase)), channel_length/(m(ibase) - 4*s_final(ibase)), 1000);
        pdf_TOF{ibase} = (channel_length^1.5./(x{ibase}.^2*s(ibase)*(2*pi*dx(ibase))^0.5)).*exp(-channel_length*(channel_length./x{ibase} - m(ibase)).^2/(2*s(ibase)^2*dx(ibase)));
    end
    
    % plot time of flight probability distribution functions
    plot(x{1}/1000, pdf_TOF{1}, "g", "linewidth", 1,x{2}/1000, pdf_TOF{2},"b","linewidth", 1, x{3}/1000, pdf_TOF{3},"m","linewidth", 1, x{4}/1000, pdf_TOF{4},"r","linewidth", 1);
    legend('dCMP','dTMP','dGMP','dAMP');
    xlabel(['Time of flight over ' num2str(channel_length/1000) ' um (us)'])
    ylabel('Time of flight PDF');
    print -dpng -color time_flight.png
    %save TOF distribution data
    for ibase = 1:nbases
        outfile = ['TOFdist_' num2str(sepeff) '_' base{ibase} '.dat'];
        outdata = [x{ibase} pdf_TOF{ibase}];
        save(outfile,'outdata','-ascii');
    end
    
    %~sepdist = sep_dist_calc_norm(vel{1}(:,1),vel{2}(:,1),vel{3}(:,1),vel{4}(:,1),nsubblcks*nblcks,blcklen_short,sepeff);
    
 end

