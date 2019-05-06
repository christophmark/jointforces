def linear(stiffness):
    return {'K_0': stiffness, 'D_0': 1e30, 'L_S': 1e30, 'D_S': 1e30}


def custom(K_0, D_0, L_S, D_S):
    return {'K_0': K_0, 'D_0': D_0, 'L_S': L_S, 'D_S': D_S}


collagen06 = {'K_0': 447, 'D_0': 0.0008, 'L_S': 0.0075, 'D_S': 0.033}
collagen12 = {'K_0': 1645, 'D_0': 0.0008, 'L_S': 0.0075, 'D_S': 0.033}
collagen24 = {'K_0': 5208, 'D_0': 0.0008, 'L_S': 0.0075, 'D_S': 0.033}
