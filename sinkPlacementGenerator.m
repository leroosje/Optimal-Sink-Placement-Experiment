N = 1000;
X = rand(1,N+4)*1000;
Y = rand(1,N+4)*1000;
processednodes = zeros(1,N+4);
processednodes(1001) = 1;
processednodes(1002) = 1;
processednodes(1003) = 1;
processednodes(1004) = 1;
newprocessed = zeros(1,N+4);
newprocessed(1001) = 1;
newprocessed(1002) = 1;
newprocessed(1003) = 1;
newprocessed(1004) = 1;
R = 50;
C = 0;
Iter = 1;

neighbor = zeros(N+4,N+4);
for m = 1:N+4
    neighbor(m,:) = ((X - X(m)).^2 + (Y-Y(m)).^2 < R*R);
    neighbor(m,m) = 0;
end

costs = double.empty; %Resets the matrix of costs so output is clear
limiter = 15; %The limiter is the dimension of the grid formed on the graph
for s = 2:limiter %Runs for each s by s grid
sD = 1000/s; %Actual dimension of each sector
ss = s*s; %Number of sectors total
P = 10; %Penalty for disconnected nodes
%Used to bound the sectors
xmin = 0;
xmax = xmin + sD;
ymin = 0;
ymax = ymin + sD;
density = zeros(s.^2);
%Populate the density of a given sector by iterating through them
for k = 0:s-1
for i = 1:s
for j = 1:1000
if(xmin<=X(j) && X(j)<=xmax && ymin<=Y(j) && Y(j)<=ymax)
density(k*s + i) = density(k*s + i) + 1;
end
end
xmin = xmin+sD;
xmax = xmax+sD;
end
ymin = ymin+sD;
ymax = ymax+sD;
xmin = 0;
xmax = xmin + sD;
end
%Rank the sectors in terms of density by marking them in a two dimensional
%matrix(Rank, Sector#)
rank = zeros(ss,ss);
curr = 1;
for i = 1:ss
for j = 1:ss
if(density(j) > density(curr))
curr = j;
end
end
rank(i,curr) = 1;
density(curr) = -1;
end
%Climbs through the ranked sectors
for L = 0:(ss) - 4
for i = 0:3
for j = 1:ss
if(rank((ss)-i-L,j)~=0)
k = j;
end
end
%Places the sinks in the proper location
tempX = sD * (mod(k-1,s)) + sD/2;
tempi = floor((k-1)/s);
tempY = sD * tempi + sD/2;
X(N + i + 1) = tempX;
Y(N + i + 1) = tempY;
end
%Find the neighbors of the newly placed sinks
for m = 1001:1004
neighbor(m,:) = ((X-X(m)).^2 + (Y-Y(m)).^2<R*R);
neighbor(m,m) = 0;
end
%Prepare for running the cost
Iter = 1;
C = 0;
processednodes = zeros(1,1004)
processednodes(1001) = 1;
processednodes(1002) = 1;
processednodes(1003) = 1;
processednodes(1004) = 1;
newprocessed = zeros(1,1004);
newprocessed(1001) = 1;
newprocessed(1002) = 1;
newprocessed(1003) = 1;
newprocessed(1004) = 1;
%figure; %Output for testing
while(sum(newprocessed)>0)
preprocessed = processednodes;
newlyprocessed = double.empty;
for j = 1:N+4
if(newprocessed(j) == 1)
newlyprocessed = [newlyprocessed j];
end
end
for k = newlyprocessed
neighbors = double.empty;
for i = 1:N+4
if(processednodes(i) == 0 && neighbor(k,i) == 1)
neighbors = [neighbors i];
end
end
for j = 1:length(neighbors)
n = neighbors(j);
C = C + Iter; %Iter represents the number hop that the network is forming on, so for each node connected at this hop, it is added to cost
%plot([X(k) X(n)],[Y(k) Y(n)],'color',[0,0,0]) %Output for testing
%hold on %Output for testing
processednodes(n) = 1;
end
end
Iter = Iter + 1;
newprocessed = processednodes - preprocessed;
end
for i = 1:N+4
if(processednodes(i) == 0)
C = C + P; %if a node wasn't processed, it is isolated. Add penalty to cost.
end
end
costs = [costs C]; %The total cost is held onto in this matrix
%plot(X(1001),Y(1001),'color',[1,0,0],'marker','.'); hold on; %Output
%plot(X(1002),Y(1002),'color',[1,0,0],'marker','.'); hold on; %Output
%plot(X(1003),Y(1003),'color',[1,0,0],'marker','.'); hold on; %Output
%plot(X(1004),Y(1004),'color',[1,0,0],'marker','.'); hold on; %Output
end
costs = [costs -1]; %Seperate limiter runs with a -1 for readability
end