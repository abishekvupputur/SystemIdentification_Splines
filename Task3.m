EKF;
%% Step 1
Ints=1000;
alpha=Z_k1k1(1,Ints:N);
beta=Z_k1k1(2,Ints:N);
x=[alpha;beta]';
y=Cm(Ints:N);
%% Step 2
V=[-0.02529 0.1483;
   -0.1242 0.1792;
   -0.1236 0.1229;];
tri=delaunayTriangulation(V);
[bc]=bsplinen_cart2bary(V,x);
index=(bc(:,1)>0 & bc(:,2)>0 & bc(:,3)>0);
bc = bc.*(bc(:,1)>0 & bc(:,2)>0 & bc(:,3)>0);
bc = bc(any(bc,2),:);
xp = bsplinen_bary2cart(V,bc);
Y=y.*index;
Y= Y(any(Y,2));
triplot(tri)
hold on
scatter(xp(:,1),xp(:,2));
d=2;
for i=1:size(bc,1)
    l=1;
for k0=d:-1:0
    for k1=d:-1:0
        for k2=d:-1:0
            k=k0+k1+k2;
            if k==d
                [k0 k1 k2];
                B(i,l)=factorial(d)/(factorial(k0)*factorial(k1)*factorial(k2))*bc(i,1)^k0*bc(i,2)^k1*bc(i,3)^k2;
                l=l+1;
            end
        end
    end
end
end
%% LSE
c=inv((B'*B))*B'*Y;
p=B*c;
plot(p)
hold on
plot(Y)
hold off
figure;
tri=delaunay(xp);
tr = triangulation(tri,xp(:,1),xp(:,2),p);
trisurf(tr,'EdgeColor', 'none')
hold on
plot3(xp(:,1),xp(:,2),Y, '.k')
