%% Main
function dy = X_dot(t,X)

disp(t)
global G density Face n

R = X(1:3);
V = X(4:6);


%% Loop for force of each tetahedron
af = 0;
for f = 1:n
    
    Facef = Face(f);
    
    nf = Facef.nf;
    nij = Facef.nij;
    njk = Facef.njk;
    nki = Facef.nki;
    
    Ri = Facef.Pi - R';
    Rj = Facef.Pj - R';
    Rk = Facef.Pk - R';
    ri = norm(Ri);
    rj = norm(Rj);
    rk = norm(Rk);
    
    wf = 2*atan(dot(Ri,cross(Rj,Rk))/(ri*rj*rk + ...
        ri*(dot(Rj,Rk)) + rj*(dot(Rk,Ri)) + rk*(dot(Ri,Rj))));
    
    Lij = log((ri + rj + Facef.eij)/(ri + rj - Facef.eij));
    Ljk = log((rj + rk + Facef.ejk)/(rj + rk - Facef.ejk));
    Lki = log((rk + ri + Facef.eki)/(rk + ri - Facef.eki));
    
    af = af + nf*(wf*dot(nf,Ri) -...
        Lij*dot(nij,Rj) - Ljk*dot(njk,Rj) - Lki*dot(nki,Rk));
    
end

a = G*density*af;
dy = [V;a'];
end