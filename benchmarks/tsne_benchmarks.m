addpath('/Users/george/Research_Local/FIt-SNE')
sigma_0 = .0001;
cluster_num = 10;
Ns =   [1E4 1E5 5E5 1E6];
fi_times = zeros(length(Ns), 1);
bh_times = zeros(length(Ns), 1);
dim = 5;
%% Save the Ps
for Ni=1:length(Ns)
    rng(3)
    N = Ns(Ni)
    cluster_size = N/cluster_num;
    X_clusters = [];
    for i =1:cluster_num
        X_clusters = [X_clusters; 0.001*mvnrnd(rand(dim,1)', diag(sigma_0*ones(dim,1)), cluster_size)];
    end
    % Save input similarities
    clear opts
    opts.stop_lying_iter = 1;
    opts.max_iter = 1;
    opts.rand_seed = 3;
    opts.load_affinities = 'save';
    fast_tsne(X_clusters, opts);
    % FIt-SNE
    clear opts
    opts.rand_seed = 3;
    opts.load_affinities = 'load';
    tic 
    fast_tsne(X_clusters, opts);
    fi_times(Ni) = toc
    % BH
    clear opts
    opts.rand_seed = 3;
    opts.nbody_algo = 'bh';
    opts.load_affinities = 'load';
    tic 
    fast_tsne(X_clusters, opts);
    bh_times(Ni) = toc
end 
%save('bh_fi_nbody.mat', 'bh_times','fi_times')



%% Test the input similarities


sigma_0 = .0001;
cluster_num = 10;
Ns =   [1E4 1E5 5E5 1E6];
vptree_times = zeros(length(Ns), 2);
vptreeMT_times = zeros(length(Ns), 2);
annoy_times = zeros(length(Ns), 2);
dims = [25 50 100];
%% Save the Ps
rng(3)
for Ni=1:length(Ns)
    for Di=1:length(dims)
        N = Ns(Ni)
        dim = dims(Di)
        cluster_size = N/cluster_num;
        X_clusters = [];
        for i =1:cluster_num
            X_clusters = [X_clusters; 0.001*mvnrnd(rand(dim,1)', diag(sigma_0*ones(dim,1)), cluster_size)];
        end
        % Annoy
        clear opts
        opts.stop_lying_iter = 1;
        opts.max_iter = 1;
        tic
        fast_tsne(X_clusters, opts);
        annoy_times(Ni,Di) = toc;
        %VPTree MT
        clear opts
        opts.stop_lying_iter = 1;
        opts.max_iter = 1;
        opts.knn_algo = 'vptree';
        tic
        fast_tsne(X_clusters, opts);
        vptreeMT_times(Ni,Di) = toc;
        %VPTree 1t       
        clear opts
        opts.stop_lying_iter = 1;
        opts.max_iter = 1;
        opts.knn_algo = 'vptree';
        opts.nthreads = 1;
        tic
        fast_tsne(X_clusters, opts);
        vptree_times(Ni,Di) = toc;
        [vptree_times vptreeMT_times annoy_times]
        save('vptree_annoy_is.mat', 'vptree_times', 'vptreeMT_times', 'annoy_times')
    end
end 




