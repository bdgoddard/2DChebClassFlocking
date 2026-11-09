function r = makeGridSlit(optsPhys)
% put particles in a slit [0,L].
% used to get initial configuration for sampling

nParticles=sum(optsPhys.nParticlesS);
L = optsPhys.L;

if(isfield(optsPhys,'sigma'))
    sigma = optsPhys.sigma; 
    sigma = sigma(1);
else
    sigma=1;
end

allowedSlit = [0 + sigma; L - sigma]; % add a bit of leeway at the edges

r = linspace(allowedSlit(1), allowedSlit(2), nParticles);

r = r(:);

end

