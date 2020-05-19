function [c,ceq] = CED_Constraint(x,Returns,Target,PathLength, i, Alpha)

%CED Constraint
Portfolio = x*Returns(i-2*PathLength:i, :)';

ced = CED_Faster(Portfolio', PathLength, Alpha);

c = ced - Target;

ceq = [];

end


