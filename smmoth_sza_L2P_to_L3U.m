dL=0.02;
p =[-dL   90+dL/2]; 
q =[ dL -180-dL/2]; 
N=180/dL;
M=360/dL;


filename_L3U='20171007203000-STAR-L3U_GHRSST-SSTsubskin-VIIRS_NPP-ACSPO_V2.50B04-v02.0-fv01.0.nc';
filename_L2='20171007221000-STAR-L2P_GHRSST-SSTsubskin-VIIRS_NPP-ACSPO_V2.50B04-v02.0-fv01.0.nc';

QL=double(ncread(filename_L3U,'quality_level'))';
sza_L2=double(ncread(filename_L2,'satellite_zenith_angle')');

[I,J]=size(sza_L2);
S_q=zeros(I,J);
for i=1:I
    d=diff(sza_L2(i,:));
    ind=find(d~=0);
    mid_ind=[1 round((ind(1:end-1)+ind(2:end))/2) J];
    S_q(i,:)=interp1(mid_ind,sza_L2(i,mid_ind),[mid_ind(1):mid_ind(end)]);
end

sst_dtime=double(ncread(filename_L2,'sst_dtime')');
T_UTC=str2num(filename_L2(9:10))+str2num(filename_L2(11:12))/60+sst_dtime/(60*60);

lat_L2=double(ncread(filename_L2,'lat')');
lon_L2=double(ncread(filename_L2,'lon')');
ii=round((lat_L2-p(2))/p(1));
jj=round((lon_L2-q(2))/q(1));

H=zeros(N,M); T=zeros(N,M); S=zeros(N,M);
for i=1:I
    for j=1:J
        H(ii(i,j),jj(i,j))=H(ii(i,j),jj(i,j))+1;
        T(ii(i,j),jj(i,j))=T(ii(i,j),jj(i,j))+T_UTC(i,j);
        S(ii(i,j),jj(i,j))=S(ii(i,j),jj(i,j))+S_q(i,j);
    end
end

w=5;
f=fspecial('gaussian',2*w+1,w/2);
sza_L3_smooth=NaN*ones(N,M); time_L3_smooth=NaN*ones(N,M);
for n=1:N
    for m=1:M
        if  isfinite(QL1(n,m))==1
            h_win=H(max(1,n-w):min(N,n+w),max(1,m-w):min(M,m+w));
            s_win=S(max(1,n-w):min(N,n+w),max(1,m-w):min(M,m+w));
            t_win=T(max(1,n-w):min(N,n+w),max(1,m-w):min(M,m+w));
            ind=find(h_win~=0);
            if length(ind)>=1
                if length(h_win(:))==(2*w+1)^2
                    sza_L3_smooth(n,m)=sum(f(ind).*(s_win(ind)./h_win(ind)))/sum(f(ind));
                    time_L3_smooth(n,m)=sum(f(ind).*(t_win(ind)./h_win(ind)))/sum(f(ind));
                else
                    sza_L3_smooth(n,m)=mean(s_win(ind)./h_win(ind));
                    time_L3_smooth(n,m)=mean(t_win(ind)./h_win(ind));
                end
            end
        end
    end
end

% save sza_L3_smooth as sza
% save time_L3_smooth sa time