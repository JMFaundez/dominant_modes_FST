function U = MeanFlow_BL(y,Uin)
% U = MeanFlow_readBF(y,Uin,yin)
% This function interpolate the inpit flow filed on Chebychev points
% 
% y: array containing coordinates of grid in normal direction
% Uin: input meanflow field
% yin: wall-normal coordinate of the input meanflow
%
% Output is the structre array 'U' which contains 
% U(i).u, U(i).v, U(i).w: streamwise, normal and spanwise velocity
% U(i).uy, U(i).vy, U(i).wy: normal derivatives of u, v and w 
% U(i).ux, U(i).wx: streawwise derivatives of u and w 
%
% (c) Ardeshir Hanifi & David Tempelmann, 2014
%

% set global reference scales
method = 'pchip';

% rescale meanflow such that all quantities are based on reference scales of
% the station x=1 and interpolate onto Gauss-Lobatto Grid

Nstation=length(Uin);
U=struct();
for i = 1:Nstation
    yin=Uin(i).y;
    U(i).u=interp1(yin,Uin(i).u,y,method,Uin(i).u(end));
    U(i).w=interp1(yin,Uin(i).w,y,method,Uin(i).w(end));
    U(i).v=interp1(yin,Uin(i).v,y,method,Uin(i).v(end));
 
    U(i).uy=interp1(yin,Uin(i).uy,y,method,0);
    U(i).wy=interp1(yin,Uin(i).wy,y,method,0);
    U(i).vy=interp1(yin,Uin(i).vy,y,method,Uin(i).vy(end));

    U(i).ux=interp1(yin,Uin(i).ux,y,method,Uin(i).ux(end));
    U(i).wx=interp1(yin,Uin(i).wx,y,method,Uin(i).wx(end));
    U(i).vx=interp1(yin,Uin(i).vx,y,method,Uin(i).vx(end));  
    
    for k=1:length(y)
        if y(k)>yin(end)
            U(i).v(k)=Uin(i).v(end)+Uin(i).vy(end)*(y(k)-yin(end));
        end
    end
end


end
