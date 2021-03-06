a=10;
b = 10;
d = a*b;
nsamples = 2000;

% projection operator for constraining weights to 2D lattice
non_zero = eye(d);
for i = 1:a
    for j = 1:b
        idx1 = (i-1)*b+j;
        if i==1 || i==a || j==1 || j==b
            if i==1
                non_zero(idx1, i*b+j) = 1;
                non_zero(idx1, (a-1)*b+j) = 1;
                if j>1
                    non_zero(idx1, (i-1)*b+j-1) = 1;
                end
                if j<b
                    non_zero(idx1, (i-1)*b+j+1) = 1;
                end
            end
            if i==a
                non_zero(idx1, j) = 1;
                non_zero(idx1, (i-2)*b+j) = 1;
                if j>1
                    non_zero(idx1, (i-1)*b+j-1) = 1;
                end
                if j<b
                    non_zero(idx1, (i-1)*b+j+1) = 1;
                end

            end
            if j==1
                non_zero(idx1, (i-1)*b+j+1) = 1;
                non_zero(idx1, (i-1)*b+b) = 1;
                if i>1
                    non_zero(idx1, (i-2)*b+j) = 1;
                end
                if i<a
                    non_zero(idx1, i*b+j) = 1;
                end
            end
            if j==b
                non_zero(idx1, (i-1)*b+j-1) = 1;
                non_zero(idx1, (i-1)*b+1) = 1;
                if i>1
                    non_zero(idx1, (i-2)*b+j) = 1;
                end
                if i<a
                    non_zero(idx1, i*b+j) = 1;
                end

            end
        else
            non_zero(idx1, (i-2)*b+j) = 1;
            non_zero(idx1, i*b+j) = 1;
            non_zero(idx1, (i-1)*b+j-1) = 1;
            non_zero(idx1, (i-1)*b+j+1) = 1;
        end
    end
end

weights = randn(d)/sqrt(4) * 3.0;
weights = (weights + weights')/2;
weights = weights.*non_zero;
weights = weights - diag(diag(weights));
biases = -sum( weights );
biases = biases(:);

state_init = (rand(d,nsamples) < 0.5);
%state_init(:) = 0;
state = state_init;

% GO!
iters = 10;
stateholder = {};
for i=1:iters
    tic();
    [state, cprobs] = sw_allall_bm(state, weights, biases, 50);
    stateholder = [stateholder, state];
    toc()
end

[state_thorough, cprobs] = sw_allall_bm(state, weights, biases, 500);
Clast = corrcoef( state_thorough' );
fprintf('\n\n\n\n');
for i = 1:iters
    figure(i);
    subplot( 1, 2, 1 );
    C = corrcoef( stateholder{i}' );
    imagesc( C); colorbar;
    title( sprintf( 'C step %d', i ) );
    subplot( 1, 2, 2 );
    imagesc( C-Clast ); colorbar;
    title( sprintf( 'C-C_{true} step %d', i ) );    
    fprintf( 'step %d L2 error %f\n', i, sum((C(:) - Clast(:)).^2) );
end