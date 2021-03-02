//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "FieldMappingDiffusionMaterial.h"

registerMooseObject("ExampleApp", FieldMappingDiffusionMaterial);

template <>
InputParameters
validParams<FieldMappingDiffusionMaterial>()
{
  InputParameters params = validParams<Material>();

  // UserObjectName is the MOOSE type used for getting the name of a UserObject from the input file
  params.addRequiredParam<UserObjectName>(
      "field_mapping_userobject",
      "The name of the UserObject that is going to be computing the "
      "average value of a variable on each block");

  return params;
}

FieldMappingDiffusionMaterial::FieldMappingDiffusionMaterial(const InputParameters & parameters)
  : Material(parameters),

    // Declare that this material is going to provide a Real
    // valued property named "diffusivity" that Kernels can use.
    _diffusivity(declareProperty<Real>("diffusivity")),

    // When getting a UserObject from the input file pass the name
    // of the UserObjectName _parameter_
    // Note that getUserObject returns a _const reference_ of the type in < >
    _field_mapping_value(getUserObject<FieldMappingValue>("field_mapping_userobject"))
{
}

void
FieldMappingDiffusionMaterial::computeQpProperties()
{
  // We will compute the diffusivity based on the average value of the variable on each block.

  // We'll get that value from a UserObject that is computing it for us.

  // To get the current block number we're going to query the "subdomain_id()" of the current
  // element

  const std::vector<Real> & vector_test = _field_mapping_value.getValue();

  _diffusivity[_qp] = 0.5*vector_test[0];
}
