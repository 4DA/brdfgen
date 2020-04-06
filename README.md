# BRDF Generator
This is a simple program that generates Smith GGX BRDF lookup tables for the split sum approximation of the UE4-based PBR pipeline. 

This is rewrite in Julia of https://github.com/HectorMF/BRDFGenerator

# Bootstrapping julia installation

# Usage

# Algorithm
```
For each pixel (x, y):
   Compute roughness (0 to 1.0) based on pixel x coordinate.
   Compute NoV (0 to 1.0) based on pixel y coordinate.
   Set view as float3(sqrt(1.0 - NoV * NoV), 0, NoV).
   Set normal as float3(0, 0, 1). 
   For each sample:
        Compute a Hammersely coordinate.
        Integrate number of importance samples for (roughness and NoV).
        Compute reflection vector L
        Compute NoL (normal dot light)
        Compute NoH (normal dot half)
        Compute VoH (view dot half)
        
        If NoL > 0
          Compute the geometry term for the BRDF given roughness squared, NoV, NoL
          Compute the visibility term given G, VoH, NoH, NoV, NoL
          Compute the fresnel term given VoH.
          Sum the result given fresnel, geoemtry, visibility.
   Average result over number of samples.
```
