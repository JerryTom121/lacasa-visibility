%% build mex files if necessary
mex conway.c
mex visibility.c

%% dimensions
len = 5000;

%% make Conway series
tic
y = conway(len);
toc

%% build visibility matrix
tic
v = visibility(y);
toc

%% plot series & log-log degree distribution
subplot(121)
plot(y);
subplot(122)
[n, x] = hist(sum(v), 100);
loglog(x, n);