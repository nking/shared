# some methods are from keras code  (Apache 2.0 License https://github.com/keras-team/keras/blob/v2.14.0/LICENSE)
# https://github.com/keras-team/keras/blob/v2.14.0/keras/utils/dataset_utils.py
# https://github.com/keras-team/keras/blob/v2.14.0/keras/utils/image_dataset.py
#
# else author: Nichole King

import tempfile
import numpy as np
import time
import os
from PIL import Image
from tempfile import TemporaryDirectory
import tensorflow as tf
from tensorflow import keras
import glob
import random
import shutil

# from keras code  (Apache 2.0 License https://github.com/keras-team/keras/blob/v2.14.0/LICENSE)
# https://github.com/keras-team/keras/blob/v2.14.0/keras/utils/dataset_utils.py
def paths_and_labels_to_dataset(image_paths, image_size, num_channels, labels,
    label_mode, num_classes, interpolation, crop_to_aspect_ratio=False,):
    """Constructs a dataset of images and labels."""
    # TODO(fchollet): consider making num_parallel_calls settable
    path_ds = tf.data.Dataset.from_tensor_slices(image_paths)
    args = (image_size, num_channels, interpolation, crop_to_aspect_ratio)
    img_ds = path_ds.map(
        lambda x: load_image(x, *args), num_parallel_calls=tf.data.AUTOTUNE
    )
    img_ds = img_ds.reshape(tf.reshape(s))

    if label_mode:
        label_ds = labels_to_dataset(
            labels, label_mode, num_classes
        )
        img_ds = tf.data.Dataset.zip((img_ds, label_ds))
    return img_ds

# from keras code  (Apache 2.0 License https://github.com/keras-team/keras/blob/v2.14.0/LICENSE)
# https://github.com/keras-team/keras/blob/v2.14.0/keras/utils/dataset_utils.py
def labels_to_dataset(labels, label_mode, num_classes):
    """Create a tf.data.Dataset from the list/tuple of labels.

    Args:
      labels: list/tuple of labels to be converted into a tf.data.Dataset.
      label_mode: String describing the encoding of `labels`. Options are:
      - 'binary' indicates that the labels (there can be only 2) are encoded as
        `float32` scalars with values 0 or 1 (e.g. for `binary_crossentropy`).
      - 'categorical' means that the labels are mapped into a categorical
        vector.  (e.g. for `categorical_crossentropy` loss).
      num_classes: number of classes of labels.

    Returns:
      A `Dataset` instance.
    """
    label_ds = tf.data.Dataset.from_tensor_slices(labels)
    if label_mode == "binary":
        label_ds = label_ds.map(
            lambda x: tf.expand_dims(tf.cast(x, "float32"), axis=-1),
            num_parallel_calls=tf.data.AUTOTUNE,
        )
    elif label_mode == "categorical":
        label_ds = label_ds.map(
            lambda x: tf.one_hot(x, num_classes),
            num_parallel_calls=tf.data.AUTOTUNE,
        )
    return label_ds

## from keras code  (Apache 2.0 License https://github.com/keras-team/keras/blob/v2.14.0/LICENSE)
## https://github.com/keras-team/keras/blob/v2.14.0/keras/utils/image_dataset.py
def load_image(
    path, image_size, num_channels, interpolation, crop_to_aspect_ratio=False
):
    """Load an image from a path and resize it."""
    img = tf.io.read_file(path)
    img = tf.image.decode_image(
        img, channels=num_channels, expand_animations=False
    )
    if crop_to_aspect_ratio:
        img = tf.keras.preprocessing.image.smart_resize(
            img, image_size, interpolation=interpolation
        )
    else:
        img = tf.image.resize(img, image_size, method=interpolation)
    img.set_shape((image_size[0], image_size[1], num_channels))
    return img

def get_subset_train_val_lists_0(data_dir: str, n_train: int, fraction_val: float, seed : int, verbose:int=0) :
    '''
    given the path to the top of the specific dataset directory
    of directory substructure:
       e.g. give path to cifar10 for this:
           cifar10
              train
                 airplanes
                 cars
                 ...
               test
                  airplanes
                  cars
                  ...
    and given the number of train samples wanted, the method extracts
    all train filepaths and selects n_train from them randomly.
    From the remaining file_paths, the code selects the validation
    list randomly.
    Then returns as a tuple (train_files, train_labels, val_files, val_labels, class_names)

    Note that other dataset directory substructures need a different method
    to be created for them.
    '''

    if fraction_val > 1 or fraction_val < 0:
        raise ValueError("fraction_val must be >= 0 and <= 1")

    _ = set(glob.glob(data_dir + '/train/**/*'))

    class_names = {p.split('/')[-1]: i for i, p in enumerate(sorted(glob.glob(data_dir + '/train/**')))}

    n_val = int(n_train*0.8)

    train_files = random.sample(sorted(_), n_train)
    _.difference_update(set(train_files))
    val_files = random.sample(sorted(_), n_val)

    if verbose > 0:
        print(f'n_class_names = {len(class_names)}')
        print(f'n_train file paths extracted = {len(train_files)}')
        print(f'n_val file paths extracted = {len(val_files)}')

    train_labels = np.zeros((n_train,), dtype="int32")
    for i, fl in enumerate(train_files):
        # './tempcifar10/train/airplanes/img0015.png'
        # directory above filename.  parse string for last 2 directory separators in os indep manner
        directories = fl.split("/")
        label = class_names[directories[-2]]
        if verbose > 1:
            print(f'filename={fl}\nlabel={label}  {directories[-2]}')
        train_labels[i] = label

    val_labels = np.zeros((n_val,), dtype="int32")
    for i, fl in enumerate(val_files):
        directories = fl.split("/")
        label = class_names[directories[-2]]
        if verbose > 1:
            print(f'filename={fl}\nlabel={label}  {directories[-2]}')
        val_labels[i] = label

    return (train_files, train_labels, val_files, val_labels, class_names)

# NOTE: DO NOT USE.  dataset internal tensor shapes are not correct for the model
# ValueError: Input 0 of layer "model" is incompatible with the layer: expected shape=(None, 32, 32, 3),
# found shape=(32, 32, 3)
def get_subset_train_val_datasets_0(image_size: tuple, num_channels: int, crop_to_aspect_ratio: bool, batch_size: int, data_dir:str, n_train:int, fraction_val: float, seed:int, verbose:int=0):
    '''
    given the path to the top of the specific dataset directory
    of directory substructure:
       e.g. give path to cifar10 for this:
           cifar10
              train
                 airplanes
                 cars
                 ...
               test
                  airplanes
                  cars
                  ...
    and given the number of train samples wanted, the method extracts
    all train filepaths and selects n_train from them randomly.
    From the remaining file_paths, the code selects the validation
    list randomly.
    The class labels are extracted from the file paths.
    Builds tf.data.Dataset structures from the file paths and labels.

    Note that other dataset directory substructures need a different method
    to be created for them.
    
    returns [train_dataset, val_dataset, class_names]
    '''

    train_files, train_labels, val_files, val_labels, class_names = get_subset_train_val_lists_0(data_dir, n_train,
                                                          fraction_val=fraction_val, seed=seed, verbose=verbose)
    train_dataset = paths_and_labels_to_dataset(
        image_paths=train_files,
        image_size=image_size,
        num_channels=num_channels,
        labels=train_labels,
        label_mode='int',
        num_classes=len(class_names),
        interpolation='bilinear',
        crop_to_aspect_ratio=crop_to_aspect_ratio,
    )

    val_dataset = paths_and_labels_to_dataset(
        image_paths=val_files,
        image_size=image_size,
        num_channels=num_channels,
        labels=val_labels,
        label_mode='int',
        num_classes=len(class_names),
        interpolation='bilinear',
        crop_to_aspect_ratio=crop_to_aspect_ratio,
    )

    train_dataset.class_names = class_names
    val_dataset.class_names = class_names
    train_dataset.file_paths = train_files
    val_dataset.file_paths = val_files

    return [train_dataset, val_dataset, class_names]

def get_subset_test_list_0(data_dir: str, n_test: int, seed: int, verbose:int = 0):
    '''
    Args:
        data_dir - top level directory of project.  e.g. ./cifar10
        n_test - number of test files.  if None, all tests are used
    '''
    _ = set(glob.glob(data_dir + '/test/**/*'))

    class_names = {p.split('/')[-1]: i for i, p in enumerate(sorted(glob.glob(data_dir + '/test/**')))}

    if n_test is None:
        n_test = len(_)
        test_files = _
    else:
        test_files = random.sample(sorted(_), n_test)

    if verbose > 0:
        print(f'n_class_names = {len(class_names)}')
        print(f'n_test file paths extracted = {len(test_files)}')

    test_labels = np.zeros((n_test,), dtype="int32")
    for i, fl in enumerate(test_files):
        # './tempcifar10/test/airplanes/img0015.png'
        directories = fl.split("/")
        label = class_names[directories[-2]]
        if verbose > 1:
            print(f'filename={fl}\nlabel={label}  {directories[-2]}')
        test_labels[i] = label

    return (test_files, test_labels, class_names)

def copy_subset_to_TMP(data_dir: str, n_train: int, n_test: int, fraction_val: float, seed: int, verbose: int = 0):
    '''
    given the path to the top of the specific dataset directory
    of directory substructure:
       e.g. give path to cifar10 for this:
           cifar10
              train
                 airplanes
                 cars
                 ...
               test
                  airplanes
                  cars
                  ...
    and given the number of train samples wanted, the method extracts
    all train filepaths and selects n_train from them randomly.
    From the remaining file_paths, the code selects the validation
    list randomly.
    if n_test is None, uses all test files.

    random samples all files to copy n_train and n_train*fraction_val to a new temporary
    directory to be loaded using keras methods
    '''

    train_files, train_labels, val_files, val_labels, class_names = get_subset_train_val_lists_0(
        data_dir, n_train, fraction_val, seed, verbose)

    test_files, test_labels, test_class_names = get_subset_test_list_0(data_dir, n_test, seed, verbose)

    #persist to file system
    tmp_dir = data_dir + '_TMP'
    n_del = len(data_dir)

    if not os.path.exists(tmp_dir):
        for fl in train_files:
            fl2 = tmp_dir + fl[n_del:]
            dir_path = fl2[0:fl2.rindex('/')]
            os.makedirs(dir_path, exist_ok=True)
            shutil.copyfile(fl, fl2)

        for fl in val_files:
            fl2 = tmp_dir + fl[n_del:]
            dir_path = fl2[0:fl2.rindex('/')]
            os.makedirs(dir_path, exist_ok=True)
            shutil.copyfile(fl, fl2)

        for fl in test_files:
            fl2 = tmp_dir + fl[n_del:]
            dir_path = fl2[0:fl2.rindex('/')]
            os.makedirs(dir_path, exist_ok=True)
            shutil.copyfile(fl, fl2)

    return tmp_dir

