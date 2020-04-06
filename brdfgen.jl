using Images
using Colors
using FileIO
using LinearAlgebra
using ImageView
using FileIO
using Printf

const Vec2 = Array{Float64, 1}
const Vec3 = Array{Float64, 1}

# function RadicalInverse_VdC(bits::Int16) 
# {
#     bits = (bits << 16u) | (bits >> 16u);
#     bits = ((bits & 0x55555555u) << 1u) | ((bits & 0xAAAAAAAAu) >> 1u);
#     bits = ((bits & 0x33333333u) << 2u) | ((bits & 0xCCCCCCCCu) >> 2u);
#     bits = ((bits & 0x0F0F0F0Fu) << 4u) | ((bits & 0xF0F0F0F0u) >> 4u);
#     bits = ((bits & 0x00FF00FFu) << 8u) | ((bits & 0xFF00FF00u) >> 8u);
#     return float(bits) * 2.3283064365386963e-10;
# }

function VanDerCorpus(n::Int64, base::Int64)::Float64
    invBase = 1.0 / float(base);
    denom   = 1.0;
    result  = 0.0;

    for i = 0:31
        if(n > 0)
            denom   = mod(float(n), 2.0);
            result += denom * invBase;
            invBase = invBase / 2.0;
            n       = convert(Int64, floor(float(n) / 2.0));
        end
    end

    return result;
end

function Hammersley(i::Int64, N::Int64)::Vec2
    return [convert(Float64,i) / convert(Float64, N), VanDerCorpus(i, 2)]
end

function ImportanceSampleGGX(Xi::Vec2, roughness::Float64, N::Vec3)::Vec3
    a = roughness*roughness;

    phi = 2.0 * pi * Xi[1];
    cosTheta = sqrt((1.0 - Xi[2]) / (1.0 + (a*a - 1.0) * Xi[2]));
    sinTheta = sqrt(1.0 - cosTheta*cosTheta);

    #  from spherical coordinates to cartesian coordinates
    H = [cos(phi) * sinTheta, sin(phi) * sinTheta, cosTheta]
         
    #  from tangent-space vector to world-space sample vector
    up = abs(N[3]) < 0.999 ? [0.0, 0.0, 1.0] : [1.0, 0.0, 0.0];
    tangent = normalize(cross(up, N));
    bitangent = cross(N, tangent);

    sampleVec = tangent * H[1] + bitangent * H[2] + N * H[3];
    return normalize(sampleVec);
end

function GeometrySchlickGGX(NdotV::Float64, roughness::Float64)::Float64
    a = roughness;
    k = (a * a) / 2.0;

    nom = NdotV;
    denom = NdotV * (1.0 - k) + k;

    return nom / denom;
end


function GeometrySmith(roughness::Float64, NoV::Float64, NoL::Float64)::Float64
    ggx2 = GeometrySchlickGGX(NoV, roughness);
    ggx1 = GeometrySchlickGGX(NoL, roughness);
    return ggx1 * ggx2;
end


function integrateBRDF(NdotV::Float64, roughness::Float64, samples::Int64)::Vec2
    V = [sqrt(1.0 - NdotV * NdotV), 0.0, NdotV]
    A = 0.0;
    B = 0.0;
    N = [0.0, 0.0, 1.0];

    for i = 0:samples-1
	Xi = Hammersley(i, samples);
	H = ImportanceSampleGGX(Xi, roughness, N);
	L = normalize(2.0 * dot(V, H) * H - V);

	NoL = max(L[3], 0.0);
	NoH = max(H[3], 0.0);
	VoH = max(dot(V, H), 0.0);
	NoV = max(dot(N, V), 0.0);

	if NoL > 0.0
	    G = GeometrySmith(roughness, NoV, NoL);

	    G_Vis = (G * VoH) / (NoH * NoV);
	    Fc = (1.0 - VoH)^5.0;

	    A += (1.0 - Fc) * G_Vis;
	    B += Fc * G_Vis;
	end
    end

    return [A / convert(Float64, samples), B / convert(Float64, samples)]
end

function main()
    samples = 1000
    size = 512
    progress::Int64 = 0

    data = rand(RGB{Float32}, size,size)
    for x = 0:size-1
        for y = 0:size-1
	    NoV = (y + 0.5) * (1.0 / size);
	    roughness = (x + 0.5) * (1.0 / size);
            rg = integrateBRDF(NoV, roughness, samples);
            data[x+1,y+1] = RGB{Float32}(rg[1], rg[2], 0.0)
        end

        progressNew = convert(Int64, ceil(100 * x / size))

        if (progressNew > progress)
            progress = progressNew
            if (progress % 10 == 0)
                @printf("%d%% ", progress)
            end
        end
    end

    @printf("\n")

    
    imshow(data)
    # to check with reference pictures
    rgb24 = convert.(RGB24, data)
    save(File(format"PNG", "brdfLUT.png"), colorview(RGB,rgb24))
end

main()
