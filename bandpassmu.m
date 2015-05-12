function y=bandpassmu(x, sf, lco, hco)
% bandpass filter
% lco: low cut-off frequency (Hz)
% hco: high cut-off frequency (Hz)
% order: order of filter
% sf: sampling frequency (Hz)

order=6;

n=round(order/2);
wn=[lco hco]/(sf/2);
[b,a]=butter(n, wn);
y=filtfilt(b,a,x);

