#ifndef __itkGradientMagnitudeImageToPointSetMetric_h
#define __itkGradientMagnitudeImageToPointSetMetric_h

#include "itkImageToPointSetMetric.h"
#include "itkCovariantVector.h"
#include "itkPoint.h"
#include <itkSignedMaurerDistanceMapImageFilter.h>


namespace itk
{
/** \class GradientMagnitudeImageToPointSetMetric
 * \brief Computes similarity between pixel values of a point set and
 * intensity values of an image.
 *
 * This metric computes the correlation between point values in the fixed
 * point-set and pixel values in the moving image. The correlation is
 * normalized by the autocorrelation values of both the point-set and the
 * moving image. The spatial correspondence between the point-set and the image
 * is established through a Transform. Pixel values are taken from the fixed
 * point-set. Their positions are mapped to the moving image and result in
 * general in non-grid position on it.  Values at these non-grid position of
 * the moving image are interpolated using a user-selected Interpolator.
 *
 * \ingroup RegistrationMetrics
 * \ingroup ITKRegistrationCommon
 */
template< class TFixedImage, class TMovingPointSet >
class ITK_EXPORT MinimumValueImageToPointSetMetric:
  public ImageToPointSetMetric< TFixedImage, TMovingPointSet >
{
public:
  typedef MinimumValueImageToPointSetMetric                Self;
  typedef ImageToPointSetMetric< TFixedImage, TMovingPointSet > Superclass;
  typedef SmartPointer< Self >                                  Pointer;
  typedef SmartPointer< const Self >                            ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(MinimumValueImageToPointSetMetric,
               ImageToPointSetMetric);

  /** Types transferred from the base class */
  typedef typename Superclass::RealType                RealType;
  typedef typename Superclass::TransformType           TransformType;
  typedef typename Superclass::TransformPointer        TransformPointer;
  typedef typename Superclass::TransformParametersType TransformParametersType;
  typedef typename Superclass::TransformJacobianType   TransformJacobianType;
  typedef typename Superclass::GradientPixelType       GradientPixelType;

  typedef typename Superclass::MeasureType                MeasureType;
  typedef typename Superclass::DerivativeType             DerivativeType;
  typedef typename Superclass::MovingPointSetType         MovingPointSetType;
  typedef typename Superclass::FixedImageType             FixedImageType;
  typedef typename Superclass::MovingPointSetConstPointer MovingPointSetConstPointer;
  typedef typename Superclass::FixedImageConstPointer     FixedImageConstPointer;

  typedef typename Superclass::PointIterator     PointIterator;
  typedef typename Superclass::PointDataIterator PointDataIterator;
  typedef typename Superclass::InputPointType    InputPointType;
  typedef typename Superclass::OutputPointType   OutputPointType;

	typedef Image<float, itkGetStaticConstMacro(FixedImageDimension)>  DistanceMapType;
	typedef itk::SmartPointer<DistanceMapType>     DistanceMapPointer;

  /** Get the derivatives of the match measure. */
  void GetDerivative(const TransformParametersType & parameters,
                     DerivativeType & Derivative) const;

  /**  Get the value for single valued optimizers. */
  MeasureType GetValue(const TransformParametersType & parameters) const;

  /**  Get value and derivatives for multiple valued optimizers. */
  void GetValueAndDerivative(const TransformParametersType & parameters,
                             MeasureType & Value, DerivativeType & Derivative) const;

  /** Set/Get SubtractMean boolean. If true, the sample mean is subtracted
   * from the sample values in the cross-correlation formula and
   * typically results in narrower valleys in the cost fucntion.
   * Default value is false. */
  itkSetMacro(SubtractMean, bool);
  itkGetConstReferenceMacro(SubtractMean, bool);
  itkBooleanMacro(SubtractMean);

	void Initialize(void);
protected:
  MinimumValueImageToPointSetMetric();
  virtual ~MinimumValueImageToPointSetMetric() {}
  void PrintSelf(std::ostream & os, Indent indent) const;

private:
  MinimumValueImageToPointSetMetric(const Self &); //purposely not
                                                            // implemented
  void operator=(const Self &);                             //purposely not
                                                            // implemented

  bool m_SubtractMean;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkMinimumValueImageToPointSetMetric.hxx"
#endif

#endif
