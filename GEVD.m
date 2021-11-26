function [sig_yy,sig_nn,sig_ss,d,Q,Qh,V] = GEVD(Rnn,Ryy)
    % Computing the MWF filter using a generalised eigenvalue
    % decomposition of the correlation matrices.
    [V,d] = eig(Ryy,Rnn,'vector');
    [d,ind] = sort(d,'descend','ComparisonMethod','real');
    V = V(:,ind);
    Qh = pinv(V);
    sig_yy = V'*Ryy *V;
    sig_yy = [sig_yy(1,1) 0;0 sig_yy(2,2)];
    sig_nn = V'*Rnn*V;
    sig_nn = [sig_nn(1,1) 0;0 sig_nn(2,2)];
    Q= Qh';
    sig_ss = sig_yy(1,1) - sig_nn(1,1) ; 
    sig_ss = [sig_ss 0; 0 0];
end