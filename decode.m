function digits = decode(signal)
% Q9:
% implement touchtone decoder using STFT

win = [zeros(384,1); ones(512,1); zeros(640,1)];  % this window is tailor made for the frequencies that are relevant to touchtone signal.
hopSize = 1536;  % assume Fs = 4096, this hop size will give us exactly one tone for each time slot.
F = 128;  % the number of frequencies for the STFT. this value was chosen ampiricaly.

S = STFT(signal, win, hopSize, F, 1);  % calculating the STFT
S_size = size(S);

Fs = 4096;
F_vec = linspace(-Fs/2, Fs/2, S_size(1));  % calculating the frequencies vector.

% calculating the index in the frequencies vector for each one of the
% frequencies that arte relevant to the touchtone signal.
[~,ind_1209]=min(abs(F_vec-1209));
[~,ind_1336]=min(abs(F_vec-1336));
[~,ind_1477]=min(abs(F_vec-1477));
[~,ind_697]=min(abs(F_vec-697));
[~,ind_770]=min(abs(F_vec-770));
[~,ind_852]=min(abs(F_vec-852));

% For dealing with additive gaussian noise, we'll make an average of the STFT at
% the frequency range of each relevant frequency (for example, for 1209Hz take
% average of 1177Hz - 1241Hz)
S_high = abs(S([ind_1209, ind_1336, ind_1477],:))+abs(S([ind_1209-1, ind_1336-1, ind_1477-1],:))+abs(S([ind_1209+1, ind_1336+1, ind_1477+1],:));
S_low = abs(S([ind_697, ind_770, ind_852],:))+abs(S([ind_697-1, ind_770-1, ind_852-1],:))+abs(S([ind_697+1, ind_770+1, ind_852+1],:));
% average
S_high = S_high/3;
S_low = S_low/3;

% calculating the high frequency and the low frequency of each tone by
% choosing the frequecy with the highest power in each time slot
[~,I_high] = max(abs(S_high),[],1);
[~,I_low] =  max(abs(S_low),[],1);

% if we look at the table matching between frequecies and digits, we can 
% define I_high as the culomns and I_low as the rows
% (column 1 is 1209Hz and row 1 is 697Hz)
% we can see in the table that : digit = column_num + row_num*3
digits = I_high + (I_low-1)*3;

end