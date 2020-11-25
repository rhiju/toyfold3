function C_eff = C_eff_gaussian_chain_func( D, gaussian_variance,  force_combined_gaussian_approximation  )
% C_eff = C_eff_gaussian_chain_func( D, gaussian_variance,  force_combined_gaussian_approximation )
%
% Compute cost of chain closure of chain of rods with lengths D and
%   Gaussian chains with total gaussian variance gaussian_variance.
%
% Adapted from core/scoring/func/GaussianChainFunc.cc in Rosetta,
%  developed by R. Das, 2016.
%
% INPUTS
%  D = number or array of numbers with rod lengths (Angstroms).
%  gaussian_variance = total variance in Gaussian chain connectors (Angstrom^2)
%                        Rosetta loop_close assumes 5 Å^2 for each
%                        nucleotide linkage.
% Optional input:
%  force_combined_gaussian_approximation
%            = Force using combined Gaussian approximation even for 1,2,3,4
%                  rod examples.
% OUTPUT
%  C_eff = effective concentration for closing chain (M)
%
% (C) R. Das, Stanford University 2020



if ~exist( 'force_combined_gaussian_approximation','var'); force_combined_gaussian_approximation = 0; end

score = 0;
if length(D) > 4 | force_combined_gaussian_approximation
    score = gaussian_chain_general_func( D, gaussian_variance );
else
    if length( D ) == 1
        score = gaussian_chain_single_func( D, gaussian_variance );
    elseif length( D ) == 2
        score = gaussian_chain_double_func( D, gaussian_variance );
    elseif length( D ) == 3
        score = gaussian_chain_triple_func( D, gaussian_variance );
    else
        assert( length( D ) == 4 );
        score = gaussian_chain_quadruple_func( D, gaussian_variance );
    end
end

C_eff = exp(-score)/(6.022e23/1e27);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rosetta NOTES...
% ---------------------------------
% Closure energies for loop cycles.
% ---------------------------------
%
%          D1
%       ------~
%      ~       | D4
%     ~        |
%    ~         ~
%    |        ~
% D2 |       /
%    |      /  D3
%     ~~~~~/
%
%  D1, D2, D3, D4 mark rigid joints, with lengths D1, ... D4.
%  Squiggles (~~~) mark gaussian chains,
%   with total variance of gaussian_variance.
%
%  In simplest case (one joint D, gaussian variance sigma^2 ):
%
%  P( closure ) = (capture volume) x ( 1 / 2 pi sigma^2 )^(3/2)
%                                  x exp( - D / 2 sigma^2 )
%
%  [For folks used to thinking about mean end-to-end distance L,
%        sigma^2 = L^2 / 3 ].
%
%  Following includes explicit formulae for one, two, three, and four joints.
%  There is also a general formula which would be straightforward to implement,
%  which I have worked out (and is basically implemented in the four-joint
%  function GaussianChainQuadrupleFunc.cc)
%
% See core/scoring/loop_graph/LoopGraph.cc for use case.
%
% More notes at:
%   https:%docs.google.com/a/stanford.edu/file/d/0B6gpwdY_Bgd0OHdzVWJGTHBvTzg/edit
%
% Rhiju, Sep. 2013
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function score = gaussian_chain_single_func( z, gaussian_variance );
score = 1.5 * log( 2 * pi * gaussian_variance );
loop_harmonic = (z * z)/ ( 2 * gaussian_variance );
score = score + loop_harmonic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function score = gaussian_chain_double_func( D, gaussian_variance );
D1 = D(1);
D2 = D(2);

score = 1.5 * log( 2*pi );
score = score +	log( 2 );
score = score + 0.5 * log( gaussian_variance );

term1 = exp( -( D1 - D2) * (D1 - D2)/ (2 * gaussian_variance ) );
term2 = exp( -( D1 + D2) * (D1 + D2)/ (2 * gaussian_variance ) );
loop_energy = -log( ( term1 - term2 ) / (D1 * D2 ));
score = score + loop_energy;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function score = gaussian_chain_triple_func( D, gaussian_variance );
D1 = D(1);
D2 = D(2);
D3 = D(3);

score = log(16*pi);

s = sqrt( 2 * gaussian_variance );
term0 = erf( ( D1 + D2 + D3 )/ s );
term1 = erf( (-D1 + D2 + D3 )/ s );
term2 = erf( ( D1 - D2 + D3 )/ s );
term3 = erf( ( D1 + D2 - D3 )/ s );
loop_energy = - log( ( -term0 + term1 + term2 + term3 ) / (D1 * D2 * D3));
score = score + loop_energy;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function score = gaussian_chain_quadruple_func( D, gaussian_variance );

D1 = D(1);
D2 = D(2);
D3 = D(3);
D4 = D(4);

score = log( 32.0 );
score = score + 1.5 * log( pi );
score = score -0.5 * log( 2 * gaussian_variance );

s = sqrt( 2 * gaussian_variance );
L_sum = 0.0;
for sgn2 = [-1,1]
    for sgn3 = [-1,1]
        for sgn4 = [-1,1]
            D_sum = D1 + ( sgn2 * D2 ) + ( sgn3 * D3 ) + ( sgn4 * D4 );
            D_sum = D_sum/s;
            sgn = -1 * sgn2 * sgn3 * sgn4;
            L_sum = L_sum +  sgn * L( D_sum );
        end
    end
end

loop_energy = - log( L_sum / (D1 * D2 * D3 * D4) );

score = score + loop_energy;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = L(x)
% note that this is the integral of erf( x ).
% up to a sqrt( pi ) factor.
% that's the key to a general expression.
val = (sqrt( pi ) * x * erf( x )) + exp( -1.0 * x * x);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function score = gaussian_chain_general_func( D, gaussian_variance )

N = length(D);% number of fixed-distance segments.

% this is a 'prefactor' in the probability.
score = ( N + 2 ) * log( 2 );
score = score + log( pi );
score = score + ( N - 3 ) * -(1/2) * log( 2 * gaussian_variance );

% Need to get all combinations of signs, e.g., if N = 4, we need the 16 multiplets, (+1,+1,+1,+1), (-1,+1,+1,+1) ...
num_terms = 2^N;
s = sqrt( 2 * gaussian_variance );
L_sum = 0;
for k = 1:num_terms    
    sgn = zeros(1,N);
    sgn_prod = 1;
    D_sum = 0.0;
    count = k;
    for i = 1:N
        sgn( i ) = 2 * mod(count,2)  - 1; % converts 0,1 to -1,1
        count = floor(count/ 2); % next digit, in binary
        D_sum = D_sum + sgn(i) * D(i);
        sgn_prod = sgn_prod * -sgn( i );
    end
    D_sum = D_sum/s;
    L_sum = L_sum + (-1) * sgn_prod * erfc_integral( D_sum, N - 3 );
end
        
loop_energy = -log( L_sum );
for i = 1:N
    loop_energy = loop_energy + log( D(i) );
end

score = score + loop_energy;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = erfc_integral(x, N)
assert( N >= -3 );
val = NaN;
if ( N == -3 ) 
    % -3th integral (i.e., third derivative ) of erfc
    val = ( 2.0 / sqrt( pi ) ) * 2 * ( -1 + 2 * x * x ) * exp( -1.0 * x * x );
elseif ( N == -2 ) 
    % -2th integral (i.e., second derivative ) of erfc
    val = ( 2.0 / sqrt( pi ) ) * 2 * x * exp( -1.0 * x * x );
elseif ( N == -1 ) 
    % -1th integral (i.e., first derivative ) of erfc
    val = ( 2.0 / sqrt( pi ) ) * exp( -1.0 * x * x );
elseif ( N == 0 ) 
    % zeroth integral of erfc, i.e. erfc
    val = erfc( x );
else
    % N-th integral can be derived from a recurrence if we know N-1-th and N-2-th integrals
    % (Note: easy to derive via integration by parts.)
    % Abramowitz and Stegun, p. 299
    val = ( ( -x / N ) * erfc_integral( x, N - 1 ) + ( 1 / ( 2.0 * N ) ) * erfc_integral( x, N - 2 ) );
end





