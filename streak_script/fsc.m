function dy = fsc(t,y,mfac)
% dy = fsc(t,y,mfac) sets up the functions coresponding to
% Falkner-Skan-Cooke simililarity soloution
%
% Input:
%   t: variable required by ode45
%   y: array containing function values corresponding to differential equation for similarity solution
%   mfac: Falkne-Skan velocity parameter U=x^mfac 
% 
% Output:
%   dy: derivative of function y
%
% (c) Ardeshir Hanifi & David Tempelmann, 2014
%


dy=zeros(5,1);
betaH=2*mfac/(mfac+1);
dy(1)=y(2);
dy(2)=y(3);
dy(2)=y(3);
dy(3)=-(y(1)*y(3)+betaH*(1-y(2)*y(2)));
dy(4)=y(5);
dy(5)=-y(1)*y(5);
end

