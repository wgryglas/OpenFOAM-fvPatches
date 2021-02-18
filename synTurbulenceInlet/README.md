# Examples of usage of this BC inside 0/U file

## Uniform mean velocity and turbulence parameters

```
INLET_AB
{
    type	synTurbulenceInlet; // esensiall - selects this BC

    nu 1.8e-5; //kinematic viscosity

    stats true; //optional, false by defult. Decides if some fields statistics should be printed

    referenceType fixed; //optional, fixed by default. This controls how referenceField is defined. In this case fixed mean that the "refrenceField" will not change

    referenceField	uniform (4.68403 0.004750831 0.0); //or you can use "referenceField $internalField" for convinience. This controls the mean veolocity

    dxmin 0.005; // optional, controls the minimal mesh size used for max wave number calculation. If not present BC will visit all edges in mesh and choos the shortes one size.

    properties
    {
        type fixed;   //optional, fixed by default. Fixed means that all turbulence parameters are the same across BC and time.
        turbIntensity 0.065; //required, turbulent intensity - variable used for Urms calculations from refrenceField
        turbScale 0.02; //required, used for other turb. prop. computation

        timeScale 3e-3; // optional, if not specified will be calculated using k and epsilon: k / (eps + SMALL)
    }

    //value $internalField; // optional, will be used for computing fluctations when solver will continue calculations from given time step.
}
```


## Interpolated mean velocity and fixed turbulence parameters

In case of interpolated values you need to provide data inside following case directory:
```
/constant/boundaryData/<boundary name>/
         --> points - file with list of points in the form:
             58
             (
             (0.06 0.00005 -0.0225)
             (0.06 0.00010 -0.0225)
             (0.06 0.00015 -0.0225)
             ...
             )
         --> 0/U  - file with list of velocity values in the form:
             58
             (
             (0.06 0.00005 -0.0225)
             (0.06 0.00010 -0.0225)
             (0.06 0.00015 -0.0225)
             ...
             )
```
Then in the BC spec. you don't need to define "referenceField" entry.
Note, the root dictionary (content of { ... } ) is the root of definition for MappedFile utility for "referenceField" interpolation.
You can place there more properties related to the MappedFile utility:
\table
    Property     | Description                            | Required | Default
    setAverage   | Use average value                      | no       | false
    perturb      | Perturb points for regular geometries  | no       | 1e-5
    points       | Name of points file                    | no       | points
    fieldTable   | Alternative field name to sample       | no       | this field name
    mapMethod    | Type of mapping                        | no       | planarInterpolation
    offset       | Offset to mapped values                | no       | Zero
\endtable


```
INLET_AB
{
    type	synTurbulenceInlet; // esensiall - selects this BC

    nu 1.8e-5; //kinematic viscosity

    stats true; //optional, false by defult. Decides if some fields statistics should be printed

    referenceType interpolated; //now required because fixed is default. This controls how referenceField is defined. In this case "refrenceField" will be interpolated basing on boundaryData

    dxmin 0.005; // optional, controls the minimal mesh size used for max wave number calculation. If not present BC will visit all edges in mesh and choos the shortes one size.

    properties
    {
        type fixed;   //optional, fixed by default. Fixed means that all turbulence parameters are the same across BC and time.
        turbIntensity 0.065; //required, turbulent intensity - variable used for Urms calculations from refrenceField
        turbScale 0.02; //required, used for other turb. prop. computation

        timeScale 3e-3; // optional, if not specified will be calculated using k and epsilon: k / (eps + SMALL)
    }

    //value $internalField; // optional, will be used for computing fluctations when solver will continue calculations from given time step.
}
```

## Interpolated mean velocity and interpolate turbulence parameters

In case of the all variables interpolated you need to provide data inside following case directory:
```
/constant/boundaryData/<boundary name>/
         --> points - file with list of points in the form:
             58
             (
             (0.06 0.00005 -0.0225)
             (0.06 0.00010 -0.0225)
             (0.06 0.00015 -0.0225)
             ...
             )
         --> 0/U  - file with list of velocity values in the form:
             58
             (
             (0.06 0.00005 -0.0225)
             (0.06 0.00010 -0.0225)
             (0.06 0.00015 -0.0225)
             ...
             )
         -->0/k - file with list of point values for kinetic turbulence in the form:
             58
             (
             0.06
             0.06
             0.06
             ...
             )
         -->0/omega - file with list of point values for the specific rate of dissipation in the form:
             58
             (
             0.06
             0.06
             0.06
             ...
             )
```
Then in the BC spec. you don't need to define "referenceField" entry. Also you need to place 2 new dictionaries under
properties:
 k {}
 omega {}
The root dictionary (content of { ... } ) is the root of definition for MappedFile utility for "referenceField" interpolation.
While "k{}" and "omega{}" is dictionary for MappedFile utility for both k,omega respectively.

Inside dictionaries you can place there more properties related to the MappedFile utility:
\table
    Property     | Description                            | Required | Default
    setAverage   | Use average value                      | no       | false
    perturb      | Perturb points for regular geometries  | no       | 1e-5
    points       | Name of points file                    | no       | points
    fieldTable   | Alternative field name to sample       | no       | this field name
    mapMethod    | Type of mapping                        | no       | planarInterpolation
    offset       | Offset to mapped values                | no       | Zero
\endtable

If you wil leave dictionaries empty you will operate on defaults (as above).

```
INLET_AB
{
    type	synTurbulenceInlet; // esensiall - selects this BC

    nu 1.8e-5; //kinematic viscosity

    stats true; //optional, false by defult. Decides if some fields statistics should be printed

    referenceType interpolated; //now required because fixed is default. This controls how referenceField is defined. In this case "refrenceField" will be interpolated basing on boundaryData

    dxmin 0.005; // optional, controls the minimal mesh size used for max wave number calculation. If not present BC will visit all edges in mesh and choos the shortes one size.

    properties
    {
       type interpolated;

       k  //required, you might leave {...} content empty, defualts for MappedFile will be used for k variable
       {
       }

       omega //required, you might leave {...} content empty, defualts for MappedFile will be used for omega variable
       {
       }

       timeScale 3e-3; // optional, if not specified will be calculated using k and epsilon: k / (eps + SMALL)
       //or by interpolation table, similarly to k and omega
       //timeScale
       //{
       //}
    }
}
```

## Fixed mean velocity and interpolated turbulence parameters

Last option, just a mixture of previous setups:

```
INLET_AB
{
    type	synTurbulenceInlet; // esensiall - selects this BC

    nu 1.8e-5; //kinematic viscosity

    stats true; //optional, false by defult. Decides if some fields statistics should be printed

    referenceType fixed; //optional, fixed by default. This controls how referenceField is defined. In this case fixed mean that the "refrenceField" will not change

    referenceField	uniform (4.68403 0.004750831 0.0); //or you can use "referenceField $internalField" for convinience. This controls the mean veolocity

    dxmin 0.005; // optional, controls the minimal mesh size used for max wave number calculation. If not present BC will visit all edges in mesh and choos the shortes one size.

    properties
    {
       type interpolated;

       k  //required, you might leave {...} content empty, defualts for MappedFile will be used for k variable
       {
       }

       omega //required, you might leave {...} content empty, defualts for MappedFile will be used for omega variable
       {
       }

       timeScale 3e-3; // optional, if not specified will be calculated using k and epsilon: k / (eps + SMALL)
       //or by interpolation table, similarly to k and omega
       //timeScale
       //{
       //}
    }
}
```
