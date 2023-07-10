using LinearAlgebra

q = -1;
m = 1;
a = 1;
b = 1;
ω = 0.25;
dt = 0.00001;
t = 1;
N = Int64(1e6);

x₀ = [0,0,0];
v₀ = [0,1,0];
B = [0,0,10];
Ev = [1,1,1];
E = zeros(3);
x₂ = ones(3,N);
vₙ = zeros(3,N);

for i ∈ 1:3
    E[i] = Ev[i] ⋅ (a + (b ⋅ sin(ω ⋅ dt)));
end

K = @. (q⋅E) / (2⋅m);
tᵥ = @. (B⋅q⋅dt) / (2⋅m);
sᵥ = @. (2⋅tᵥ) / (1+abs(tᵥ)^2);

vₙ[:,begin] = @. v₀;

@. x₂ *= x₀;
svₙ = vₙ;
x₁ = x₀;

@time for i ∈ 1:(N-1)
    v⁻ = @. vₙ[:,i] + (K⋅dt);
    vxt = v⁻ × tᵥ;
    v₁ = @. v⁻ + vxt;
    vxs = v₁ × sᵥ;
    v⁺ = @. v⁻ + vxs;
    vₙ[:,i+1] = @. v⁺ + (K⋅dt);
    x₂[:,i] = @. x₁ + (v₁⋅dt);
    x₁ = x₂[:,i];
end

