//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Material.h"
#include "FieldMappingValue.h"

// Forward Declarations
class FieldMappingDiffusionMaterial;

template <>
InputParameters validParams<FieldMappingDiffusionMaterial>();

class FieldMappingDiffusionMaterial : public Material
{
public:
  FieldMappingDiffusionMaterial(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

private:
  /**
   * This is the member reference that will hold the computed values
   * for the Real value property in this class.
   */
  MaterialProperty<Real> & _diffusivity;

  /**
   * A member reference that will hold onto a UserObject
   * of type FieldMappingValue for us to be able to query
   * the value of a variable on each element.
   *
   * NOTE: UserObject references are _const_!
   */
  const FieldMappingValue & _field_mapping_value;
};
