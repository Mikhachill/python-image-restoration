# coding: utf-8
import Image, numpy
from functools import partial
from math import exp


def H(k, u, v, mu=0, mv=0):
#    return 3/4
    return exp(-k * ((u - mu) ** 2 + (v - mv) ** 2) ** (5/6))


def average_smoothing(pic, width, height):
    """Осредняющее сглаживание"""
    for i in xrange(width):
        for j in xrange(height):
            average_value, local_vals_counter = 0, 0
            for local_i in xrange(i - 1, i + 2):
                for local_j in xrange(j - 1, j + 2):
                    if local_i in xrange(width) and local_j in xrange(height) and not (i == local_i and j == local_j):
                        average_value += pic[local_i, local_j]
                        local_vals_counter += 1
            average_value /= local_vals_counter
            pic[i, j] = average_value
    return pic

def inverse_filtering(img, width, height):
    """Инверсная фильтрация"""
    k = 0.002
    eps = 1E-5
    pic_a = numpy.asarray(img)
    F = numpy.fft.fftshift(numpy.fft.fft2(pic_a))
    N = 5 + 3 * numpy.random.randn(len(F), len(F[0]))
    H_values = [[H(k, i, j) for j in xrange(len(F[0]))] for i in xrange(len(F))]
    N = numpy.fft.fftshift(numpy.fft.fft2(N))
    G = [[hij * fij + nij for hij, fij, nij in zip(hi, fi, ni)] for hi, fi, ni in zip(H_values,F, N)]
    degraded_image = abs(numpy.fft.ifft2(numpy.fft.ifftshift(G)))
    H_reversed = [[1 / Hij if Hij > eps else 0 for Hij in Hi] for Hi in H_values]
#    H_reversed = []
#    for i in xrange(len(H_values)):
#        H_reversed.append([])
#        for j in xrange(len(H_values[0])):
#            if H_values[i][j] > eps:
#                H_reversed[i].append(1 / H_values[i][j])
#            else:
#                print 'zero!'
#                H_reversed[i].append(0)
#    print H_reversed[0]
#    print G[0]
    restored_image = [[gij * hrij for gij, hrij in zip(gi, hri)] for gi, hri in zip(G, H_reversed)]
    restored_image = abs(numpy.fft.ifft2(numpy.fft.ifftshift(restored_image)))
    bad_pic = Image.fromarray(degraded_image.astype(numpy.uint8))
    bad_pic.save('images/bmstu_inv_noised.png')
    good_pic = Image.fromarray(numpy.asarray(restored_image).astype(numpy.uint8))
    good_pic.save('images/bmstu_inv_restored.png')


def winear_filtering(img, width, height):
    """Винеровская фильтрация"""
    k = 0.005
    eps = 1E-5
    pic_a = numpy.asarray(img)
    F = numpy.fft.fftshift(numpy.fft.fft2(pic_a))
    noise = 5 + 3 * numpy.random.randn(len(F), len(F[0]))
    H_values = [[H(k, i, j) for j in xrange(len(F[0]))] for i in xrange(len(F))]
    N = numpy.fft.fftshift(numpy.fft.fft2(noise))
    G = [[hij * fij + nij for hij, fij, nij in zip(hi, fi, ni)] for hi, fi, ni in zip(H_values,F, N)]
    degraded_image = abs(numpy.fft.ifft2(numpy.fft.ifftshift(G)))

    Sn = [[abs(n) ** 2 for n in ni] for ni in numpy.fft.fftshift(numpy.fft.fft2(N))]
    nA = sum([sum(sni) for sni in Sn]) / (len(Sn) * len(Sn[0]))
    Sf = [[abs(n) ** 2 for n in ni] for ni in numpy.fft.fftshift(numpy.fft.fft2(pic_a))]
    fA = sum([sum(sfi) for sfi in Sf]) / (len(Sf) * len(Sf[0]))
    R = nA / fA
    R *= 0.001
    print 'R = %s' % R

    H_reversed = [[1 / Hij if Hij > eps else 0 for Hij in Hi] for Hi in H_values]
    H_squared = [[abs(Hij ** 2) for Hij in Hi] for Hi in H_values]

    W = []
    for i in xrange(len(H_values)):
        W.append([])
        for j in xrange(len(H_values[0])):
            val = H_reversed[i][j] * H_squared[i][j] / (H_values[i][j] * H_squared[i][j] + R) * G[i][j]
            W[i].append(val)

    restored_image = abs(numpy.fft.ifft2(numpy.fft.ifftshift(W)))

    bad_pic = Image.fromarray(degraded_image.astype(numpy.uint8))
    bad_pic.save('images/bmstu_winear_noised.png')
    good_pic = Image.fromarray(numpy.asarray(restored_image).astype(numpy.uint8))
    good_pic.save('images/bmstu_winear_restored.png')


if __name__ == '__main__':
    img = Image.open('images/bmstu_gray.png')
#    img = img.convert('L')
    pic = img.load()
    width, height = img.size
#    average_smoothing(pic, width, height)
    winear_filtering(img, width, height)
