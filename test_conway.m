%% build mex files if necessary
mex CFLAGS='-fPIC -std=c99 -O3 -ffast-math -march=native' conway.c
mex CFLAGS='-fPIC -std=c99 -O3 -ffast-math -march=native' visibility4.c

%% dimensions
len = 50000;

%% make Conway series
tic
y = conway(len);
toc

%% with Lacasa's method
tic
v2 = visibility2(y);
toc
%%
k2 = vis_to_spmat(v2);
figure, spy(k2)


%% with Lacasa's method, dynamic alloc
tic
[v4, ~, ~, ~] = visibility4(y);
toc
%%
k4 = vis_to_spmat(double(v4));
figure, spy(k4)

%% plot series & log-log degree distribution
subplot(121)
plot(y);
subplot(122)
[n, x] = hist(sum(k2), 100);
loglog(x, n);