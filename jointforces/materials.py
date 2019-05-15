def linear(stiffness):
    """
    Defines the stiffness for a linear elastic material.
    
    Args:
        stiffness(float) : Bulk modulus of the material in Pascal for SAENO Simulation (see [Steinwachs,2015])
    """
    return {'K_0': stiffness, 'D_0': 1e30, 'L_S': 1e30, 'D_S': 1e30}


def custom(K_0, D_0, L_S, D_S):
    """
    Defines the material properties for a custom nonlinear material.
    
    Args:
        K_0(float) : Bulk modulus of the material in Pascal for SAENO Simulation (see [Steinwachs,2015])
        D_0(float) : Buckling coefficient of the fibers for SAENO Simulation (see [Steinwachs,2015])
        L_S(float) : Onset of strain stiffening of the fibers (see [Steinwachs,2015])
        D_S(float) : Strain stiffening coefficient of the fibers (see [Steinwachs,2015])   
    """
    return {'K_0': K_0, 'D_0': D_0, 'L_S': L_S, 'D_S': D_S}


collagen06 = {'K_0': 447, 'D_0': 0.0008, 'L_S': 0.0075, 'D_S': 0.033}
collagen12 = {'K_0': 1645, 'D_0': 0.0008, 'L_S': 0.0075, 'D_S': 0.033}
collagen24 = {'K_0': 5208, 'D_0': 0.0008, 'L_S': 0.0075, 'D_S': 0.033}
