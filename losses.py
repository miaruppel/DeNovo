import numpy


def masked_spectral_distance(true, pred):
    # Note, fragment ions that cannot exists (i.e. y20 for a 7mer) must have the value  -1.
    import tensorflow.compat.v1 as tf
    import keras.backend as k
    
    epsilon = k.epsilon()
    pred_masked = ((true + 1) * pred) / (true + 1 + epsilon)
    true_masked = ((true + 1) * true) / (true + 1 + epsilon)
    pred_norm = k.l2_normalize(true_masked)
    true_norm = k.l2_normalize(pred_masked)
    product = k.sum(pred_norm * true_norm, axis=1)
    arccos = tf.acos(product)
    return 2 * arccos / numpy.pi


losses = {"masked_spectral_distance": masked_spectral_distance}


def get(loss_name):
    if loss_name in losses:
        return losses[loss_name]
    else:
        return loss_name
