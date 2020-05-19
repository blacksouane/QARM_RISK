function [c,ceq] = VolConstraint(x,Target,CovMat)

c = sqrt(252)*sqrt(x*CovMat*x') - Target;

ceq = [];

end


