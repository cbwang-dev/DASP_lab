% BEWARE!!! after running this script, the RIRs is needed to be regenerated
% for the Computer_RIRs.mat is modified as a nonn-sweetspot scenario!!!
SOE
% Answer to 1.4.10 (sweet spot) By manipulating m_pos and apply in 
% create_rirs.m
if sweetspotFlag == 1
    synth_errors=[];
    for i = 1:10
%         m_pos(:,1)=m_pos(:,1)+ones(2,1)*0.01;
        m_pos(:,2)=m_pos(:,2)+ones(2,1)*0.01;
        [RIR_sources,~]=create_rirs(m_pos,s_pos,v_pos,room_dim,rev_time,fs_RIR,Lh);
        HL=[]; 
        for j=1:J
            part_toeplitz=toeplitz(RIR_sources(:,1,j),zeros(Lg-1,1));
            HL=[HL part_toeplitz];
        end
        HR=[];  
        for j=1:J
            part_toeplitz=toeplitz(RIR_sources(:,2,j),zeros(Lg-1,1));
            HR=[HR part_toeplitz];
        end
        H = [HL;HR];
        estimated_HRTF_total=H*g;
        estimated_HRTF_left =estimated_HRTF_total(1:Lh);
        estimated_HRTF_right=estimated_HRTF_total(Lh+1:end);
        estimated_left=conv(source,estimated_HRTF_left);
        estimated_right=conv(source,estimated_HRTF_right);
        estimated_binaural_new=[estimated_left,estimated_right];
        synth_errors=[synth_errors norm(estimated_binaural-estimated_binaural_new)];
    end
    figure
    plot(1:10,synth_errors)
    title('Synthesis Error with the Movement of Both Mics')
    xlabel('cm')
    % TODO: test mic movement in another direction
end


