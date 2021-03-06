%%% this function makes the computations for CV function and draws all the necessary plots
%%% knots can be changed 

function [cvc, fdobjcv, fdobj,fdobjabscv]=cvcalc(data)


knots=1:100;
knots=knots';

y=data;

%spline object is created

r1=knots(1,1);
r2=knots(end,1);
rng = [r1,r2];
norder   =4;
nbasis   = length(knots) + norder - 2;
hgtbasis = create_bspline_basis(rng, nbasis, norder, knots);
%plot(hgtbasis); 
Lfd      = int2Lfd(2);
lambda   = 10^(-4);
hgtfdPar = fdPar(hgtbasis, Lfd, lambda);

%data is inserted inside the fdobject

[fdobj, ~, ~, ~, ~, ~, ~, ~, ~] = smooth_basis_ist(knots, y, hgtfdPar);

[stdfd]= std(fdobj);



%% CV calculation

  fd=fdobj;
  coef     = getcoef(fd);
  coefd    = size(coef);
  ndim     = length(coefd);
  if (coefd(1) == 1)
    error('Only one replication found.');
  end

  basisfd  = getbasis(fd);
  nbasis   = getnbasis(basisfd);
  rangeval = getbasisrange(basisfd);

  varbifd  = var(fd);

  %neval    = 10*nbasis + 1;
  neval    = 100;
  evalarg  = linspace(rangeval(1), rangeval(2), neval)';
  vararray = eval_bifd(varbifd, evalarg, evalarg);
  nvdim    = length(size(vararray));

  if (ndim == 2)
    stdmat  = sqrt(diag(vararray));
  else
    nvar = coefd(3);
    stdmat = zeros(neval, nvar);
    m = 0;
    for j = 1:nvar
      m = m + j;
      if (nvdim == 3)
        stdmat(:,j) = sqrt(diag(vararray(:,:,1,m)));
      else
        stdmat(:,j) = sqrt(diag(vararray(:,:,m)));
      end
    end
  end
  

meanmat=mean(y,2);
cvc=stdmat./meanmat;

%CV basis function
abscvc=abs(cvc);
[fdobjabscv, ~, ~, ~, ~, ~, ~, ~, ~] = smooth_basis_ist(knots, abscvc, hgtfdPar);

[fdobjcv, ~, ~, ~, ~, ~, ~, ~, ~] = smooth_basis_ist(knots, cvc, hgtfdPar);

%drawing plots

% fig=figure('Name', 'fdobj');
% plot_ist(fdobj);
% %saveas(fig,'Figure_fdobj.pdf');
% figure('Name', 'cv')
% plot_ist(fdobjcv);
% figure('Name', 'abscv')
% plot_ist(fdobjabscv);
% figure('Name', 'std')
% plot_ist(std(fdobj));
% axis([min(knots) max(knots) min(stdmat)-min(stdmat)*0.10 max(stdmat)+max(stdmat)*.10]) 
% figure('Name', 'mean')
% plot_ist(mean(fdobj));
% axis([min(knots) max(knots) min(meanmat)+min(meanmat)*0.10 max(meanmat)+max(meanmat)*.10]) 



% %plotting of std dev.
% 
% plot(std(fdobj1));
% hold on
% F1=plot(std(fdobj1c));
% hold on;
% F2=plot(std(fdobj1e));
% set(F1,'Color',[1 0 0]);
% set(F2,'Color',[0 1 0]);
% legend('rawdata','no outlier','N-1');
% ylim([0.35 0.60]);
% hold off
% 
% plot_ist(fdobj4);
% hold on;
% F1=plot_ist(fdobj4c);
% hold on;
% F2=plot_ist(fdobj4e);
% set(F1,'Color',[1 0 0]);
% set(F2,'Color',[0 1 0]);
% axis([min(knots) max(knots) 0 40]);
% legend('rawdata','no outlier','N-1')
% 

