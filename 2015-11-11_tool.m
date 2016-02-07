%% -*- mode:octave; indent-tabs-mode:nil -*-
%% MATH 441 project three: Make an Octave/MATLAB tool for probability
%% Author: James Starkman, jas497


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plotting and testing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = samplePMF (prob, vals)
  %% Samples a discrete PDF (sometimes called a "Probability Mass Function")
  %% Parameters:
  %%   prob: a probability vector (PMF)
  %%   vals: values associated with the PDF (the outcomes)
  %% Example: rolling a pair of dice
  %%   rolled=samplePMF([1,2,3,4,5,6,5,4,3,2,1]/36, [2,3,4,5,6,7,8,9,10,11,12]);
  i       = 1;
  cdf     = 0;
  notdone = 1;
  u       = rand(1);
  bins    = length(prob);
  %% if (length(vals) != bins || abs(sum(prob)-1) > 1e-9);
  %%   error("Bad arguments.");
  %% endif
  while (notdone && i <= bins);
    if (cdf <= u && u < cdf+prob(i));
      out = vals(i);
      notdone = 0;
    else
      cdf += prob(i++);
    endif
  endwhile
  if (notdone);
    out = vals(end);
  endif
endfunction


function out = iterativeSample (n, fxn)
  %% Calls the given sampler repeatedly.
  %% Parameters:
  %%   n:   How many samples to take
  %%   fxn: Handle to a function taking no arguments that returns a sample
  %% Returns:
  %%   The outcomes.
  %% Example: rolling 1000 pairs of dice
  %%   dist = iterativeSample(1e3, @() samplePMF([1,2,3,4,5,6,5,4,3,2,1]/36,
  %%                                               [2,3,4,5,6,7,8,9,10,11,12]));
  %%   hist(dist,11); %plot it
  if (n <= 0);
    error("Bad arguments.");
  endif
  out = zeros(n,1);
  n++;
  while (n-->1);
    out(n) = fxn();
  endwhile
endfunction


function out = hist2cdf (heights)
  %% Converts the output of hist into a CDF.  Please use more than ten bins!
  %% Example:
  %%   [heights, centers] = hist(dist, binCount);
  %%   plot(centers, hist2cdf(heights));
  bins = numel(heights);
  out  = zeros(bins, 1);
  tall = sum(heights);
  cdf  = 1;                     % doing this backwards
  heights /= tall;              % scale heights
  bins++;
  while (bins-->1);             % ensures that out(end)==1
    out(bins) = cdf;
    cdf -= heights(bins);
  endwhile                      % now, cdf < epsilon
endfunction


function out = colorpicker (iterations)
  %% Produces the "rainbow" visual effect on the plots of the distributions
  switch iterations;
    case 10;   out = "r";
    case 100;  out = "g";
    case 1000; out = "b";
    otherwise; out = "m";
  endswitch
endfunction


function histogram (heights, centers, iterations)
  %% Plots the sample output, either as colored bubbles or as a histogram
  if iterations == 1e4;
    bar(centers, heights, "hist", "facecolor", "white");     
  else
    scatter(centers, heights, colorpicker(iterations));
  endif
endfunction

function cumdist (heights, centers, iterations)
  scatter(centers, hist2cdf(heights), colorpicker(iterations));
endfunction

function outs = doHists (sampler, bins)
  %% Does the repetitive histogram work, including: scaling heights, creating a
  %% convenient output object, and doing all iterations (albeit separately).
  %% Parameters:
  %%   sampler: A closure that will be passed into `iterativeSample()'
  %%   bins:    An integer number for the second argument of `hist()'
  %% Returns:
  %%   outs:    Each column is a different number of iterations.  Each row is
  %%   structured as follows:
  %%     1: heights
  %%     2: centers
  %%     3: iterations
  outs = cell(3,4);
  for iterations = [1:length(outs)];
    its = 10^iterations;
    distrib = iterativeSample(its, sampler);
    [heights, centers] = hist(distrib, bins, 1);
    heights /= (centers(2)-centers(1));
    outs{1,iterations} = heights;
    outs{2,iterations} = centers;
    outs{3,iterations} = its;
  endfor
endfunction

function doPlots (outs, x, pdfRef, cdfRef)
  %% Does the repetitive plotting.  Uses the output object from `doHists()'.
  %% Parameters:
  %%   outs:   object from `doHists()'
  %%   x:      x-values to use for the plot
  %%   pdfRef: reference PDF to plot (will be paired with `x' for plotting).
  %%   cdfRef: reference CDF to plot (will be paired with `x' for plotting).
  figure;
  hold on;
  for sim = [outs(:,4) outs(:,1:3)];
    histogram(sim{1}, sim{2}, sim{3});
  endfor
  plot(x, pdfRef(x));
  hold off;
  
  figure;
  hold on;
  plot(x, cdfRef(x));
  for sim = outs;
    cumdist(sim{1}, sim{2}, sim{3});
  endfor
  hold off;         
endfunction


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Discrete functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = pmfBinomial (p, n)
  %% Samples the binomial distribution with n Bernoulli trials (the definition).
  %% Parameters:
  %%   p: probability of success
  %%   n: number of trials
  out = sum(rand(n,1) < p);
endfunction

function simBinomial()
  probSuccess  = 0.7;
  numberTrials = 20;
  sampler      = @() pmfBinomial(probSuccess, numberTrials);
  outs = doHists(sampler, numberTrials+1);
  
  xs = [0:1:numberTrials];
  doPlots(outs, xs,
          @(x) binopdf(x, numberTrials, probSuccess),
          @(x) binocdf(x, numberTrials, probSuccess));
endfunction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = pmfPoisson (lambda)
  %% Samples the Poisson distribution of the given parameter
  out = 0;                      % how many arrived
  pro = exp(-lambda);           % initial probability to roll under
  tot = pro;                    % total so far
  u = rand(1);                  % stopping point
  while (u > tot);
    out++;
    pro *= (lambda/out);
    tot += pro;
  endwhile
endfunction

function simPoisson()
  lambda  = 4;
  sampler = @() pmfPoisson(lambda);
  outs = doHists(sampler, 15);

  xs = [0:1:15]; %shouldn't go much higher than this for lambda=4
  doPlots(outs, xs,
          @(x) poisspdf(x, lambda),
          @(x) poisscdf(x, lambda));
endfunction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function simDice()
  PMF     = [1,2,3,4,5,6,5,4,3,2,1]/36;
  outcome = [2,3,4,5,6,7,8,9,10,11,12];
  sampler = @() samplePMF(PMF, outcome);
  outs = doHists(sampler, numel(outcome));

  xs = [2:1:12];
  analyticalCdf = hist2cdf(PMF);
  doPlots(outs, xs,
          @(x) PMF(x-1),
          @(x) analyticalCdf(x-1));
endfunction


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Continuous functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [probs,outs] = invcdfNormal (mean, standardDeviation, nSlices, maximum)
  %% norminv(x, mean, standardDeviation);
  outs  = [-maximum:2*maximum/nSlices:maximum];
  cumpr = normcdf(outs, mean, standardDeviation); % cumulative probabilities
  probs = zeros(size(outs));
  i = 1;
  probs(i) = cumpr(1);
  while (i++<length(outs));
    probs(i) = cumpr(i) - cumpr(i-1);
  endwhile
endfunction

function simNormal()
  mean    = 0;
  strdDev = 1;
  maximum = 4;
  [pr,oc] = invcdfNormal(mean, strdDev, 100, maximum);
  sampler = @() samplePMF(pr,oc); %pdfNormal(mean, strdDev);
  outs = doHists(sampler, 40);
  
  xs = [-maximum:0.01:maximum];
  doPlots(outs, xs,
          @(x) normpdf(x, mean, strdDev),
          @(x) normcdf(x, mean, strdDev));
endfunction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = pdfPareto (x, shape, minimum)
  %% Analytically computes the Pareto PDF at each element of X.  This is not an
  %% Octave/MATLAB built-in function.
  i = length(x);
  out = zeros(i,1);
  i++;
  while (i-->1);
    out(i) = shape * (minimum/x(i))^shape / x(i);
  endwhile
endfunction

function out = cdfPareto (x, shape, minimum)
  %% Analytically computes the CDF at each element of X.  This is not an
  %% Octave/MATLAB built-in function.
  i = length(x);
  out = zeros(i,1);
  i++;
  while (i-->1);
    out(i) = 1 - (minimum/x(i))^shape;
  endwhile
endfunction

function [probs,outs] = invcdfPareto (shape, minimum, nSlices, maximum)
  outs  = [minimum:(maximum-minimum)/nSlices:maximum];
  cumpr = cdfPareto(outs, shape, minimum);
  probs = zeros(size(outs));
  i = 1;
  probs(i) = cumpr(1);
  while (i++<length(outs));
    probs(i) = cumpr(i) - cumpr(i-1);
  endwhile
endfunction

function simPareto()
  param   = 3;
  minimum = 1;
  maximum = 5;
  [pr,oc] = invcdfPareto(param, minimum, 100, maximum);
  sampler = @() samplePMF(pr,oc); %pdfPareto(param, minimum);
  outs = doHists(sampler, 40);

  xs = [minimum:0.01:maximum];
  doPlots(outs, xs,
          @(x) pdfPareto(x, param, minimum),
          @(x) cdfPareto(x, param, minimum));
endfunction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [probs, outs] = invcdfGamma (k, theta, nSlices, maximum)
  %% gaminv(x, k, theta);
  outs  = [0:maximum/nSlices:maximum];
  cumpr = gamcdf(outs, k, theta);
  probs = zeros(size(outs));
  i = 1;
  probs(i) = cumpr(1);
  while (i++<length(outs));
    probs(i) = cumpr(i) - cumpr(i-1);
  endwhile
endfunction

function simGamma ()
  shape   = 2;
  scale   = 3;
  maximum = 20;
  [pr,oc] = invcdfGamma(shape, scale, 100, maximum);
  sampler = @() samplePMF(pr,oc);
  outs = doHists(sampler, 40);
  
  xs = [0:0.01:maximum];
  doPlots(outs, xs,
          @(x) gampdf(x, shape, scale),
          @(x) gamcdf(x, shape, scale));
endfunction


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  main  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function main(enable)
  %% Runs the simulations.
  %% Parameters:
  %%   enable: simple binary flags for each distribution, in the order:
  %%                                 [Binomial,Poisson,Dice,Normal,Pareto,Gamma]
  %% Warning: each selected distribution calls "figure; hold on; <plots /> hold
  %% off;" twice

  qtyDists = 6;                 % number of different distributions modelled
  if (
      length(enable) != qtyDists ||
      numel(enable) != qtyDists ||
      numel(find(enable)) + numel(find(1-enable)) != qtyDists
    );
    error("Illegal arguments.");
  else;
    if (enable(1));
      simBinomial();
    endif
    if (enable(2));
      simPoisson();
    endif
    if (enable(3));
      simDice();
    endif
    if (enable(4));
      simNormal();
    endif
    if (enable(5));
      simPareto();
    endif
    if (enable(6));
      simGamma();
    endif
  endif
endfunction

main([1,0,0,1,0,0]);

