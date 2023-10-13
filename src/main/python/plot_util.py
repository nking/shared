# adapted from
# ...

# Author: Nichole L. King

import tempfile
import numpy as np
import matplotlib.pyplot as plt
import time
import os
from PIL import Image
from tempfile import TemporaryDirectory
import tensorflow as tf
from tensorflow import keras

def plot_first_tf(ds: tf.data.Dataset, class_names: list):
    plt.figure(figsize=(6,6))
    for images, labels in ds.take(1):
        for i in range(9):
            ax = plt.subplot(3, 3, i + 1)
            plt.imshow(images[i].numpy().astype("uint8"))
            plt.title(class_names[labels[i]])

def print_ds_shape_tf(ds: tf.data.Dataset,):
    #image_batch is a tensor of the shape (32, 32, 32, 3). a batch of 32 images of shape 32x32x3
    #label_batch is a tensor of the shape (32,), these are corresponding labels to the 32 images.
    for image_batch, labels_batch in ds:
        print(f'image_batch.shape = {image_batch.shape}')
        print(f'labels_batch.shape = {labels_batch.shape}')
        break

def plot_loss_acc_tf(plot_ylabel_prefix : str, history : tf.keras.callbacks.History) :
    print(history.history.keys())
    n_keys = len(history.history.keys())
    '''
    dict_keys(['loss', 'accuracy', 'val_loss', 'val_accuracy'])
    '''
    if 'acc' in history.history and n_keys == 4:
        plot_loss_acc(plot_ylabel_prefix, history.history['loss'], history.history['acc'],
                      history.history['val_loss'], history.history['val_acc'])
    elif 'acc' in history.history and n_keys == 2:
        plot_loss_acc(plot_ylabel_prefix, history.history['loss'], history.history['acc'])
    elif 'accuracy' in history.history and n_keys == 4:
        plot_loss_acc(plot_ylabel_prefix, history.history['loss'], history.history['accuracy'],
                      history.history['val_loss'], history.history['val_accuracy'])
    elif 'accuracy' in history.history and n_keys == 2:
        plot_loss_acc(plot_ylabel_prefix, history.history['loss'], history.history['accuracy'])
    else:
        raise Exception("edit this method for the history keys")

def plot_loss_acc(plot_prefix: str, loss : list, acc : list, val_loss = None, val_acc = None) :
    '''
    plot train acc and loss and overplot val acc and loss if present.
    '''
    if (val_loss is None and val_acc is not None) or (val_loss is not None and val_acc is None):
        print(f'WARNING: val arrays are only printed if both val_loss and val_acc are not None')

    hasVal = True if (val_loss is not None and val_acc is not None) else False

    plt.figure(figsize=(6, 6))

    j = 1
    # nrows, ncols, index where index starts at 1 in the upper left corner and increases to the right
    plt.subplot(1, 2, j)

    metric_values = loss
    x_axis = [_ for _ in range(len(metric_values))]
    plt.plot(x_axis, metric_values, label='train', color='purple', marker='o', linestyle='solid', linewidth=2,
             markersize=1)
    if hasVal:
        # overplot the validation
        plt.plot(x_axis, val_loss, label='val', color='green', marker='o', linestyle='solid', linewidth=2,
                 markersize=1)
    plt.ylabel('loss')
    plt.xlabel("epoch")
    plt.legend()
    j += 1

    plt.subplot(1, 2, j)

    metric_values = acc
    x_axis = [_ for _ in range(len(metric_values))]
    plt.plot(x_axis, metric_values, label='train', color='purple', marker='o', linestyle='solid', linewidth=2,
             markersize=1)
    if hasVal:
        # overplot the validation
        plt.plot(x_axis, val_acc, label='val', color='green', marker='o', linestyle='solid', linewidth=2,
                 markersize=1)
    plt.ylabel('acc')
    plt.xlabel("epoch")
    plt.legend()
    plt.show()


def print_first_image_stats(ds):
    image_batch, labels_batch = next(iter(ds))
    first_image = image_batch[0]
    # Notice the pixel values are now in `[0,1]`.
    print(np.min(first_image), np.max(first_image))

def plot_model_activations(model : keras.Model) :
    #summary = model.summary()
    '''
    Model.summary(
        line_length=None,
        positions=None,
        print_fn=None,
        expand_nested=False,
        show_trainable=False,
        layer_range=None,
    )
    '''
    for layer in model.layers:
        #config = layer.get_config()
        #Returns the current weights of the layer, as NumPy arrays.
        weights = layer.get_weights()

