function ME = Motion_Energy_cal(pathname)

cd([pathname '/Bpod_extracted'])
% pathname = uigetdir;
% cd(pathname)

temp1=dir('Bpod_mat_*ft*');
temp2=dir('dlc_*ft*');

Bpod = readmatrix(temp1.name);
DLC  = readmatrix(temp2.name);


time = Bpod(:,1);
trial_num = Bpod(:,2);
sound = Bpod(:,4);
spout   = Bpod(:,5);

% pixel velocity %
r=3;% move mean frame
thre = 0.99;% normalize threshold by (thre*100)% percentile %

EyeL = PixelVelocity(time,DLC(:,3),DLC(:,4),r);
EyeR = PixelVelocity(time,DLC(:,1),DLC(:,2),r);
EyeRaw = EyeL + EyeR;
EyeNorm= normalizeME(EyeRaw,thre);

whisker1 = PixelVelocity(time,DLC(:,11),DLC(:,12),r);
whisker2 = PixelVelocity(time,DLC(:,13),DLC(:,14),r);
whiskerRaw = whisker1+ whisker2;
whiskerNorm= normalizeME(whiskerRaw,thre);

nose1  = PixelVelocity(time,DLC(:,5),DLC(:,6),r);
nose2  = PixelVelocity(time,DLC(:,7),DLC(:,8),r);
nose3  = PixelVelocity(time,DLC(:,9),DLC(:,10),r);
noseRaw  = nose1 + nose2 + nose3;
noseNorm = normalizeME(noseRaw,thre);

tongue1= PixelVelocity_tongue(time,DLC(:,17), DLC(:,18),r);
tongue2= PixelVelocity_tongue(time,DLC(:,19),DLC(:,20),r);
tongueRaw= tongue1 + tongue2;  
tongueNorm= normalizeME(tongueRaw,thre);


ME.time = time;
ME.trial_num = trial_num;
ME.sound = sound;
ME.spout = spout;
ME.EyeRaw = EyeRaw;
ME.whiskerRaw = whiskerRaw;
ME.noseRaw = noseRaw;
ME.tongueRaw = tongueRaw;
ME.EyeNorm = EyeNorm;
ME.whiskerNorm = whiskerNorm;
ME.noseNorm = noseNorm;
ME.tongueNorm = tongueNorm;

end

function normv = normalizeME(v,thre)

tmp= sort(v);
th = tmp(round(length(v)*thre));
normv=v./th;
normv(normv>1)=1;
end

function v = PixelVelocity(t,x,y,r)

corX = correcton(x,r);
corY = correcton(y,r);

tmpX=diff(corX);
tmpY=diff(corY);
dr=sqrt(tmpX.^2 + tmpY.^2);
dt=diff(t);

v = dr./dt*1e3;

end

function v = PixelVelocity_tongue(t,x,y,r)

corX = movmean(x,r);
corY = movmean(y,r);

nonnan_x = find(~isnan(x));
nonnan_y = find(~isnan(y));
idX = intersect(nonnan_x,find(isnan(corX)));
idY = intersect(nonnan_y,find(isnan(corY)));
corX(idX) = x(idX);
corY(idY) = y(idY);

tmpX=diff(corX);
tmpY=diff(corY);
dr=sqrt(tmpX.^2 + tmpY.^2);
dt=diff(t);

v = dr./dt*1e3;
v(isnan(v))=0;
end

function corX = correcton(x,r)

nanid = find(isnan(x));
nonnanid = setdiff(1:length(x),nanid);
nonnanVal= x(nonnanid);
corVal = zeros(length(nanid),1);
for i=1:length(nanid)
   [~,j] = min(abs(nonnanid-nanid(i)));
   corVal(i)=nonnanVal(j);
end

corX=x;
corX(nanid)=corVal;
corX = movmean(corX,r);
end