/*
 * helicoptercomp.c
 *
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * Code generation for model "helicoptercomp".
 *
 * Model version              : 11.0
 * Simulink Coder version : 9.4 (R2020b) 29-Jul-2020
 * C source code generated on : Mon Feb 28 16:45:08 2022
 *
 * Target selection: quarc_win64.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: 32-bit Generic
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#include "helicoptercomp.h"
#include "helicoptercomp_private.h"

/* Block signals (default storage) */
B_helicoptercomp_T helicoptercomp_B;

/* Continuous states */
X_helicoptercomp_T helicoptercomp_X;

/* Block states (default storage) */
DW_helicoptercomp_T helicoptercomp_DW;

/* Real-time model */
static RT_MODEL_helicoptercomp_T helicoptercomp_M_;
RT_MODEL_helicoptercomp_T *const helicoptercomp_M = &helicoptercomp_M_;
static void rate_monotonic_scheduler(void);
time_T rt_SimUpdateDiscreteEvents(
  int_T rtmNumSampTimes, void *rtmTimingData, int_T *rtmSampleHitPtr, int_T
  *rtmPerTaskSampleHits )
{
  rtmSampleHitPtr[1] = rtmStepTask(helicoptercomp_M, 1);
  rtmSampleHitPtr[2] = rtmStepTask(helicoptercomp_M, 2);
  UNUSED_PARAMETER(rtmNumSampTimes);
  UNUSED_PARAMETER(rtmTimingData);
  UNUSED_PARAMETER(rtmPerTaskSampleHits);
  return(-1);
}

/*
 *   This function updates active task flag for each subrate
 * and rate transition flags for tasks that exchange data.
 * The function assumes rate-monotonic multitasking scheduler.
 * The function must be called at model base rate so that
 * the generated code self-manages all its subrates and rate
 * transition flags.
 */
static void rate_monotonic_scheduler(void)
{
  /* Compute which subrates run during the next base time step.  Subrates
   * are an integer multiple of the base rate counter.  Therefore, the subtask
   * counter is reset when it reaches its limit (zero means run).
   */
  (helicoptercomp_M->Timing.TaskCounters.TID[2])++;
  if ((helicoptercomp_M->Timing.TaskCounters.TID[2]) > 124) {/* Sample time: [0.25s, 0.0s] */
    helicoptercomp_M->Timing.TaskCounters.TID[2] = 0;
  }
}

/*
 * This function updates continuous states using the ODE1 fixed-step
 * solver algorithm
 */
static void rt_ertODEUpdateContinuousStates(RTWSolverInfo *si )
{
  time_T tnew = rtsiGetSolverStopTime(si);
  time_T h = rtsiGetStepSize(si);
  real_T *x = rtsiGetContStates(si);
  ODE1_IntgData *id = (ODE1_IntgData *)rtsiGetSolverData(si);
  real_T *f0 = id->f[0];
  int_T i;
  int_T nXc = 4;
  rtsiSetSimTimeStep(si,MINOR_TIME_STEP);
  rtsiSetdX(si, f0);
  helicoptercomp_derivatives();
  rtsiSetT(si, tnew);
  for (i = 0; i < nXc; ++i) {
    x[i] += h * f0[i];
  }

  rtsiSetSimTimeStep(si,MAJOR_TIME_STEP);
}

/* Model output function for TID0 */
void helicoptercomp_output0(void)      /* Sample time: [0.0s, 0.0s] */
{
  /* local block i/o variables */
  real_T rtb_HILReadEncoderTimebase_o1;
  real_T rtb_HILReadEncoderTimebase_o2;
  real_T rtb_HILReadEncoderTimebase_o3;
  real_T rtb_Clock;
  real_T rtb_Backgain;
  real_T rtb_Frontgain;
  real_T rtb_Sum;
  int8_T rtAction;
  if (rtmIsMajorTimeStep(helicoptercomp_M)) {
    /* set solver stop time */
    if (!(helicoptercomp_M->Timing.clockTick0+1)) {
      rtsiSetSolverStopTime(&helicoptercomp_M->solverInfo,
                            ((helicoptercomp_M->Timing.clockTickH0 + 1) *
        helicoptercomp_M->Timing.stepSize0 * 4294967296.0));
    } else {
      rtsiSetSolverStopTime(&helicoptercomp_M->solverInfo,
                            ((helicoptercomp_M->Timing.clockTick0 + 1) *
        helicoptercomp_M->Timing.stepSize0 +
        helicoptercomp_M->Timing.clockTickH0 *
        helicoptercomp_M->Timing.stepSize0 * 4294967296.0));
    }

    {                                  /* Sample time: [0.0s, 0.0s] */
      rate_monotonic_scheduler();
    }
  }                                    /* end MajorTimeStep */

  /* Update absolute time of base rate at minor time step */
  if (rtmIsMinorTimeStep(helicoptercomp_M)) {
    helicoptercomp_M->Timing.t[0] = rtsiGetT(&helicoptercomp_M->solverInfo);
  }

  if (rtmIsMajorTimeStep(helicoptercomp_M)) {
    /* S-Function (hil_read_encoder_timebase_block): '<S4>/HIL Read Encoder Timebase' */

    /* S-Function Block: helicoptercomp/Helicopter_interface/HIL Read Encoder Timebase (hil_read_encoder_timebase_block) */
    {
      t_error result;
      result = hil_task_read_encoder
        (helicoptercomp_DW.HILReadEncoderTimebase_Task, 1,
         &helicoptercomp_DW.HILReadEncoderTimebase_Buffer[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicoptercomp_M, _rt_error_message);
      } else {
        rtb_HILReadEncoderTimebase_o1 =
          helicoptercomp_DW.HILReadEncoderTimebase_Buffer[0];
        rtb_HILReadEncoderTimebase_o2 =
          helicoptercomp_DW.HILReadEncoderTimebase_Buffer[1];
        rtb_HILReadEncoderTimebase_o3 =
          helicoptercomp_DW.HILReadEncoderTimebase_Buffer[2];
      }
    }
  }

  /* FromWorkspace: '<Root>/From Workspace' */
  {
    real_T *pDataValues = (real_T *)
      helicoptercomp_DW.FromWorkspace_PWORK.DataPtr;
    real_T *pTimeValues = (real_T *)
      helicoptercomp_DW.FromWorkspace_PWORK.TimePtr;
    int_T currTimeIndex = helicoptercomp_DW.FromWorkspace_IWORK.PrevIndex;
    real_T t = helicoptercomp_M->Timing.t[0];

    /* Get index */
    if (t <= pTimeValues[0]) {
      currTimeIndex = 0;
    } else if (t >= pTimeValues[140]) {
      currTimeIndex = 139;
    } else {
      if (t < pTimeValues[currTimeIndex]) {
        while (t < pTimeValues[currTimeIndex]) {
          currTimeIndex--;
        }
      } else {
        while (t >= pTimeValues[currTimeIndex + 1]) {
          currTimeIndex++;
        }
      }
    }

    helicoptercomp_DW.FromWorkspace_IWORK.PrevIndex = currTimeIndex;

    /* Post output */
    {
      real_T t1 = pTimeValues[currTimeIndex];
      real_T t2 = pTimeValues[currTimeIndex + 1];
      if (t1 == t2) {
        if (t < t1) {
          rtb_Clock = pDataValues[currTimeIndex];
        } else {
          rtb_Clock = pDataValues[currTimeIndex + 1];
        }
      } else {
        real_T f1 = (t2 - t) / (t2 - t1);
        real_T f2 = 1.0 - f1;
        real_T d1;
        real_T d2;
        int_T TimeIndex= currTimeIndex;
        d1 = pDataValues[TimeIndex];
        d2 = pDataValues[TimeIndex + 1];
        rtb_Clock = (real_T) rtInterpolate(d1, d2, f1, f2);
        pDataValues += 141;
      }
    }
  }

  if (rtmIsMajorTimeStep(helicoptercomp_M)) {
    /* Gain: '<S4>/Travel: Count to rad' incorporates:
     *  Gain: '<S4>/Travel_gain'
     */
    helicoptercomp_B.TravelCounttorad = helicoptercomp_P.travel_gain *
      rtb_HILReadEncoderTimebase_o1 * helicoptercomp_P.TravelCounttorad_Gain;

    /* Sum: '<Root>/Sum3' incorporates:
     *  Constant: '<Root>/deg'
     *  Gain: '<S12>/Gain'
     */
    helicoptercomp_B.Sum3 = helicoptercomp_P.Gain_Gain *
      helicoptercomp_B.TravelCounttorad + helicoptercomp_P.deg_Value;

    /* Gain: '<S4>/Pitch: Count to rad' */
    helicoptercomp_B.PitchCounttorad = helicoptercomp_P.PitchCounttorad_Gain *
      rtb_HILReadEncoderTimebase_o2;

    /* Gain: '<S9>/Gain' */
    helicoptercomp_B.Gain = helicoptercomp_P.Gain_Gain_a *
      helicoptercomp_B.PitchCounttorad;

    /* Gain: '<S4>/Elevation: Count to rad' incorporates:
     *  Gain: '<S4>/Elevation_gain'
     */
    helicoptercomp_B.ElevationCounttorad = helicoptercomp_P.elevation_gain *
      rtb_HILReadEncoderTimebase_o3 * helicoptercomp_P.ElevationCounttorad_Gain;

    /* Sum: '<Root>/Sum' incorporates:
     *  Constant: '<Root>/elavation_offset [deg]'
     *  Gain: '<S7>/Gain'
     */
    helicoptercomp_B.Sum = helicoptercomp_P.Gain_Gain_lv *
      helicoptercomp_B.ElevationCounttorad +
      helicoptercomp_P.elavation_offsetdeg_Value;
  }

  /* Clock: '<S3>/Clock' incorporates:
   *  Gain: '<S2>/Gain1'
   *  Sum: '<S5>/Sum2'
   */
  rtb_Clock -= helicoptercomp_P.Gain1_Gain * helicoptercomp_B.Gain;

  /* Sum: '<Root>/Sum1' incorporates:
   *  Constant: '<Root>/Vd_bias'
   *  Gain: '<S10>/Gain'
   *  Gain: '<S2>/Gain1'
   *  Gain: '<S5>/K_pd'
   *  Gain: '<S5>/K_pp'
   *  Sum: '<S5>/Sum3'
   *  TransferFcn: '<S4>/Pitch: Transfer Fcn'
   */
  rtb_Backgain = (helicoptercomp_P.K_pp * rtb_Clock -
                  (helicoptercomp_P.PitchTransferFcn_C *
                   helicoptercomp_X.PitchTransferFcn_CSTATE +
                   helicoptercomp_P.PitchTransferFcn_D *
                   helicoptercomp_B.PitchCounttorad) *
                  helicoptercomp_P.Gain_Gain_ae * helicoptercomp_P.Gain1_Gain *
                  helicoptercomp_P.K_pd) + helicoptercomp_P.Vd_ff;

  /* Integrator: '<S3>/Integrator' */
  /* Limited  Integrator  */
  if (helicoptercomp_X.Integrator_CSTATE >= helicoptercomp_P.Integrator_UpperSat)
  {
    helicoptercomp_X.Integrator_CSTATE = helicoptercomp_P.Integrator_UpperSat;
  } else {
    if (helicoptercomp_X.Integrator_CSTATE <=
        helicoptercomp_P.Integrator_LowerSat) {
      helicoptercomp_X.Integrator_CSTATE = helicoptercomp_P.Integrator_LowerSat;
    }
  }

  /* Sum: '<S3>/Sum' incorporates:
   *  Constant: '<Root>/elevation_ref'
   *  Gain: '<S2>/Gain1'
   */
  rtb_Sum = helicoptercomp_P.elevation_ref_Value - helicoptercomp_P.Gain1_Gain *
    helicoptercomp_B.Sum;

  /* Clock: '<S3>/Clock' incorporates:
   *  Gain: '<S2>/Gain1'
   *  Gain: '<S3>/K_ed'
   *  Gain: '<S8>/Gain'
   *  TransferFcn: '<S4>/Elevation: Transfer Fcn'
   */
  rtb_Clock = (helicoptercomp_P.ElevationTransferFcn_C *
               helicoptercomp_X.ElevationTransferFcn_CSTATE +
               helicoptercomp_P.ElevationTransferFcn_D *
               helicoptercomp_B.ElevationCounttorad) *
    helicoptercomp_P.Gain_Gain_n * helicoptercomp_P.Gain1_Gain *
    helicoptercomp_P.K_ed;

  /* Sum: '<Root>/Sum2' incorporates:
   *  Constant: '<Root>/Vs_bias'
   *  Gain: '<S3>/K_ep'
   *  Integrator: '<S3>/Integrator'
   *  Sum: '<S3>/Sum1'
   */
  rtb_Frontgain = ((helicoptercomp_P.K_ep * rtb_Sum +
                    helicoptercomp_X.Integrator_CSTATE) - rtb_Clock) +
    helicoptercomp_P.Vs_ff;

  /* Clock: '<S3>/Clock' */
  rtb_Clock = helicoptercomp_M->Timing.t[0];

  /* If: '<S3>/If' incorporates:
   *  Gain: '<S3>/K_ei'
   *  Inport: '<S6>/In1'
   */
  if (rtmIsMajorTimeStep(helicoptercomp_M)) {
    rtAction = (int8_T)!(rtb_Clock >= 2.0);
    helicoptercomp_DW.If_ActiveSubsystem = rtAction;
  } else {
    rtAction = helicoptercomp_DW.If_ActiveSubsystem;
  }

  if (rtAction == 0) {
    /* Outputs for IfAction SubSystem: '<S3>/If Action Subsystem' incorporates:
     *  ActionPort: '<S6>/Action Port'
     */
    helicoptercomp_B.In1 = helicoptercomp_P.K_ei * rtb_Sum;

    /* End of Outputs for SubSystem: '<S3>/If Action Subsystem' */
  }

  /* End of If: '<S3>/If' */

  /* Gain: '<S1>/Front gain' incorporates:
   *  Sum: '<S1>/Add'
   */
  rtb_Sum = (rtb_Backgain + rtb_Frontgain) * helicoptercomp_P.Frontgain_Gain;

  /* Saturate: '<S4>/Front motor: Saturation' */
  if (rtb_Sum > helicoptercomp_P.FrontmotorSaturation_UpperSat) {
    /* Saturate: '<S4>/Front motor: Saturation' */
    helicoptercomp_B.FrontmotorSaturation =
      helicoptercomp_P.FrontmotorSaturation_UpperSat;
  } else if (rtb_Sum < helicoptercomp_P.FrontmotorSaturation_LowerSat) {
    /* Saturate: '<S4>/Front motor: Saturation' */
    helicoptercomp_B.FrontmotorSaturation =
      helicoptercomp_P.FrontmotorSaturation_LowerSat;
  } else {
    /* Saturate: '<S4>/Front motor: Saturation' */
    helicoptercomp_B.FrontmotorSaturation = rtb_Sum;
  }

  /* End of Saturate: '<S4>/Front motor: Saturation' */

  /* Gain: '<S1>/Back gain' incorporates:
   *  Sum: '<S1>/Subtract'
   */
  rtb_Sum = (rtb_Frontgain - rtb_Backgain) * helicoptercomp_P.Backgain_Gain;

  /* Saturate: '<S4>/Back motor: Saturation' */
  if (rtb_Sum > helicoptercomp_P.BackmotorSaturation_UpperSat) {
    /* Saturate: '<S4>/Back motor: Saturation' */
    helicoptercomp_B.BackmotorSaturation =
      helicoptercomp_P.BackmotorSaturation_UpperSat;
  } else if (rtb_Sum < helicoptercomp_P.BackmotorSaturation_LowerSat) {
    /* Saturate: '<S4>/Back motor: Saturation' */
    helicoptercomp_B.BackmotorSaturation =
      helicoptercomp_P.BackmotorSaturation_LowerSat;
  } else {
    /* Saturate: '<S4>/Back motor: Saturation' */
    helicoptercomp_B.BackmotorSaturation = rtb_Sum;
  }

  /* End of Saturate: '<S4>/Back motor: Saturation' */
  if (rtmIsMajorTimeStep(helicoptercomp_M)) {
    /* S-Function (hil_write_analog_block): '<S4>/HIL Write Analog' */

    /* S-Function Block: helicoptercomp/Helicopter_interface/HIL Write Analog (hil_write_analog_block) */
    {
      t_error result;
      helicoptercomp_DW.HILWriteAnalog_Buffer[0] =
        helicoptercomp_B.FrontmotorSaturation;
      helicoptercomp_DW.HILWriteAnalog_Buffer[1] =
        helicoptercomp_B.BackmotorSaturation;
      result = hil_write_analog(helicoptercomp_DW.HILInitialize_Card,
        helicoptercomp_P.HILWriteAnalog_channels, 2,
        &helicoptercomp_DW.HILWriteAnalog_Buffer[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicoptercomp_M, _rt_error_message);
      }
    }
  }
}

/* Model update function for TID0 */
void helicoptercomp_update0(void)      /* Sample time: [0.0s, 0.0s] */
{
  if (rtmIsMajorTimeStep(helicoptercomp_M)) {
    rt_ertODEUpdateContinuousStates(&helicoptercomp_M->solverInfo);
  }

  /* Update absolute time */
  /* The "clockTick0" counts the number of times the code of this task has
   * been executed. The absolute time is the multiplication of "clockTick0"
   * and "Timing.stepSize0". Size of "clockTick0" ensures timer will not
   * overflow during the application lifespan selected.
   * Timer of this task consists of two 32 bit unsigned integers.
   * The two integers represent the low bits Timing.clockTick0 and the high bits
   * Timing.clockTickH0. When the low bit overflows to 0, the high bits increment.
   */
  if (!(++helicoptercomp_M->Timing.clockTick0)) {
    ++helicoptercomp_M->Timing.clockTickH0;
  }

  helicoptercomp_M->Timing.t[0] = rtsiGetSolverStopTime
    (&helicoptercomp_M->solverInfo);

  /* Update absolute time */
  /* The "clockTick1" counts the number of times the code of this task has
   * been executed. The absolute time is the multiplication of "clockTick1"
   * and "Timing.stepSize1". Size of "clockTick1" ensures timer will not
   * overflow during the application lifespan selected.
   * Timer of this task consists of two 32 bit unsigned integers.
   * The two integers represent the low bits Timing.clockTick1 and the high bits
   * Timing.clockTickH1. When the low bit overflows to 0, the high bits increment.
   */
  if (!(++helicoptercomp_M->Timing.clockTick1)) {
    ++helicoptercomp_M->Timing.clockTickH1;
  }

  helicoptercomp_M->Timing.t[1] = helicoptercomp_M->Timing.clockTick1 *
    helicoptercomp_M->Timing.stepSize1 + helicoptercomp_M->Timing.clockTickH1 *
    helicoptercomp_M->Timing.stepSize1 * 4294967296.0;
}

/* Derivatives for root system: '<Root>' */
void helicoptercomp_derivatives(void)
{
  XDot_helicoptercomp_T *_rtXdot;
  boolean_T lsat;
  boolean_T usat;
  _rtXdot = ((XDot_helicoptercomp_T *) helicoptercomp_M->derivs);

  /* Derivatives for TransferFcn: '<S4>/Travel: Transfer Fcn' */
  _rtXdot->TravelTransferFcn_CSTATE = 0.0;
  _rtXdot->TravelTransferFcn_CSTATE += helicoptercomp_P.TravelTransferFcn_A *
    helicoptercomp_X.TravelTransferFcn_CSTATE;
  _rtXdot->TravelTransferFcn_CSTATE += helicoptercomp_B.TravelCounttorad;

  /* Derivatives for TransferFcn: '<S4>/Pitch: Transfer Fcn' */
  _rtXdot->PitchTransferFcn_CSTATE = 0.0;
  _rtXdot->PitchTransferFcn_CSTATE += helicoptercomp_P.PitchTransferFcn_A *
    helicoptercomp_X.PitchTransferFcn_CSTATE;
  _rtXdot->PitchTransferFcn_CSTATE += helicoptercomp_B.PitchCounttorad;

  /* Derivatives for TransferFcn: '<S4>/Elevation: Transfer Fcn' */
  _rtXdot->ElevationTransferFcn_CSTATE = 0.0;
  _rtXdot->ElevationTransferFcn_CSTATE +=
    helicoptercomp_P.ElevationTransferFcn_A *
    helicoptercomp_X.ElevationTransferFcn_CSTATE;
  _rtXdot->ElevationTransferFcn_CSTATE += helicoptercomp_B.ElevationCounttorad;

  /* Derivatives for Integrator: '<S3>/Integrator' */
  lsat = (helicoptercomp_X.Integrator_CSTATE <=
          helicoptercomp_P.Integrator_LowerSat);
  usat = (helicoptercomp_X.Integrator_CSTATE >=
          helicoptercomp_P.Integrator_UpperSat);
  if (((!lsat) && (!usat)) || (lsat && (helicoptercomp_B.In1 > 0.0)) || (usat &&
       (helicoptercomp_B.In1 < 0.0))) {
    _rtXdot->Integrator_CSTATE = helicoptercomp_B.In1;
  } else {
    /* in saturation */
    _rtXdot->Integrator_CSTATE = 0.0;
  }

  /* End of Derivatives for Integrator: '<S3>/Integrator' */
}

/* Model output function for TID2 */
void helicoptercomp_output2(void)      /* Sample time: [0.25s, 0.0s] */
{
  /* (no output code required) */
}

/* Model update function for TID2 */
void helicoptercomp_update2(void)      /* Sample time: [0.25s, 0.0s] */
{
  /* Update absolute time */
  /* The "clockTick2" counts the number of times the code of this task has
   * been executed. The absolute time is the multiplication of "clockTick2"
   * and "Timing.stepSize2". Size of "clockTick2" ensures timer will not
   * overflow during the application lifespan selected.
   * Timer of this task consists of two 32 bit unsigned integers.
   * The two integers represent the low bits Timing.clockTick2 and the high bits
   * Timing.clockTickH2. When the low bit overflows to 0, the high bits increment.
   */
  if (!(++helicoptercomp_M->Timing.clockTick2)) {
    ++helicoptercomp_M->Timing.clockTickH2;
  }

  helicoptercomp_M->Timing.t[2] = helicoptercomp_M->Timing.clockTick2 *
    helicoptercomp_M->Timing.stepSize2 + helicoptercomp_M->Timing.clockTickH2 *
    helicoptercomp_M->Timing.stepSize2 * 4294967296.0;
}

/* Model output wrapper function for compatibility with a static main program */
void helicoptercomp_output(int_T tid)
{
  switch (tid) {
   case 0 :
    helicoptercomp_output0();
    break;

   case 2 :
    helicoptercomp_output2();
    break;

   default :
    break;
  }
}

/* Model update wrapper function for compatibility with a static main program */
void helicoptercomp_update(int_T tid)
{
  switch (tid) {
   case 0 :
    helicoptercomp_update0();
    break;

   case 2 :
    helicoptercomp_update2();
    break;

   default :
    break;
  }
}

/* Model initialize function */
void helicoptercomp_initialize(void)
{
  /* Start for S-Function (hil_initialize_block): '<Root>/HIL Initialize' */

  /* S-Function Block: helicoptercomp/HIL Initialize (hil_initialize_block) */
  {
    t_int result;
    t_boolean is_switching;
    result = hil_open("q8_usb", "0", &helicoptercomp_DW.HILInitialize_Card);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helicoptercomp_M, _rt_error_message);
      return;
    }

    is_switching = false;
    result = hil_set_card_specific_options(helicoptercomp_DW.HILInitialize_Card,
      "update_rate=normal;decimation=1", 32);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helicoptercomp_M, _rt_error_message);
      return;
    }

    result = hil_watchdog_clear(helicoptercomp_DW.HILInitialize_Card);
    if (result < 0 && result != -QERR_HIL_WATCHDOG_CLEAR) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helicoptercomp_M, _rt_error_message);
      return;
    }

    if ((helicoptercomp_P.HILInitialize_AIPStart && !is_switching) ||
        (helicoptercomp_P.HILInitialize_AIPEnter && is_switching)) {
      {
        int_T i1;
        real_T *dw_AIMinimums = &helicoptercomp_DW.HILInitialize_AIMinimums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AIMinimums[i1] = (helicoptercomp_P.HILInitialize_AILow);
        }
      }

      {
        int_T i1;
        real_T *dw_AIMaximums = &helicoptercomp_DW.HILInitialize_AIMaximums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AIMaximums[i1] = helicoptercomp_P.HILInitialize_AIHigh;
        }
      }

      result = hil_set_analog_input_ranges(helicoptercomp_DW.HILInitialize_Card,
        helicoptercomp_P.HILInitialize_AIChannels, 8U,
        &helicoptercomp_DW.HILInitialize_AIMinimums[0],
        &helicoptercomp_DW.HILInitialize_AIMaximums[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicoptercomp_M, _rt_error_message);
        return;
      }
    }

    if ((helicoptercomp_P.HILInitialize_AOPStart && !is_switching) ||
        (helicoptercomp_P.HILInitialize_AOPEnter && is_switching)) {
      {
        int_T i1;
        real_T *dw_AOMinimums = &helicoptercomp_DW.HILInitialize_AOMinimums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOMinimums[i1] = (helicoptercomp_P.HILInitialize_AOLow);
        }
      }

      {
        int_T i1;
        real_T *dw_AOMaximums = &helicoptercomp_DW.HILInitialize_AOMaximums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOMaximums[i1] = helicoptercomp_P.HILInitialize_AOHigh;
        }
      }

      result = hil_set_analog_output_ranges(helicoptercomp_DW.HILInitialize_Card,
        helicoptercomp_P.HILInitialize_AOChannels, 8U,
        &helicoptercomp_DW.HILInitialize_AOMinimums[0],
        &helicoptercomp_DW.HILInitialize_AOMaximums[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicoptercomp_M, _rt_error_message);
        return;
      }
    }

    if ((helicoptercomp_P.HILInitialize_AOStart && !is_switching) ||
        (helicoptercomp_P.HILInitialize_AOEnter && is_switching)) {
      {
        int_T i1;
        real_T *dw_AOVoltages = &helicoptercomp_DW.HILInitialize_AOVoltages[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOVoltages[i1] = helicoptercomp_P.HILInitialize_AOInitial;
        }
      }

      result = hil_write_analog(helicoptercomp_DW.HILInitialize_Card,
        helicoptercomp_P.HILInitialize_AOChannels, 8U,
        &helicoptercomp_DW.HILInitialize_AOVoltages[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicoptercomp_M, _rt_error_message);
        return;
      }
    }

    if (helicoptercomp_P.HILInitialize_AOReset) {
      {
        int_T i1;
        real_T *dw_AOVoltages = &helicoptercomp_DW.HILInitialize_AOVoltages[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOVoltages[i1] = helicoptercomp_P.HILInitialize_AOWatchdog;
        }
      }

      result = hil_watchdog_set_analog_expiration_state
        (helicoptercomp_DW.HILInitialize_Card,
         helicoptercomp_P.HILInitialize_AOChannels, 8U,
         &helicoptercomp_DW.HILInitialize_AOVoltages[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicoptercomp_M, _rt_error_message);
        return;
      }
    }

    if ((helicoptercomp_P.HILInitialize_EIPStart && !is_switching) ||
        (helicoptercomp_P.HILInitialize_EIPEnter && is_switching)) {
      {
        int_T i1;
        int32_T *dw_QuadratureModes =
          &helicoptercomp_DW.HILInitialize_QuadratureModes[0];
        for (i1=0; i1 < 8; i1++) {
          dw_QuadratureModes[i1] = helicoptercomp_P.HILInitialize_EIQuadrature;
        }
      }

      result = hil_set_encoder_quadrature_mode
        (helicoptercomp_DW.HILInitialize_Card,
         helicoptercomp_P.HILInitialize_EIChannels, 8U,
         (t_encoder_quadrature_mode *)
         &helicoptercomp_DW.HILInitialize_QuadratureModes[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicoptercomp_M, _rt_error_message);
        return;
      }
    }

    if ((helicoptercomp_P.HILInitialize_EIStart && !is_switching) ||
        (helicoptercomp_P.HILInitialize_EIEnter && is_switching)) {
      {
        int_T i1;
        int32_T *dw_InitialEICounts =
          &helicoptercomp_DW.HILInitialize_InitialEICounts[0];
        for (i1=0; i1 < 8; i1++) {
          dw_InitialEICounts[i1] = helicoptercomp_P.HILInitialize_EIInitial;
        }
      }

      result = hil_set_encoder_counts(helicoptercomp_DW.HILInitialize_Card,
        helicoptercomp_P.HILInitialize_EIChannels, 8U,
        &helicoptercomp_DW.HILInitialize_InitialEICounts[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicoptercomp_M, _rt_error_message);
        return;
      }
    }

    if ((helicoptercomp_P.HILInitialize_POPStart && !is_switching) ||
        (helicoptercomp_P.HILInitialize_POPEnter && is_switching)) {
      uint32_T num_duty_cycle_modes = 0;
      uint32_T num_frequency_modes = 0;

      {
        int_T i1;
        int32_T *dw_POModeValues =
          &helicoptercomp_DW.HILInitialize_POModeValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POModeValues[i1] = helicoptercomp_P.HILInitialize_POModes;
        }
      }

      result = hil_set_pwm_mode(helicoptercomp_DW.HILInitialize_Card,
        helicoptercomp_P.HILInitialize_POChannels, 8U, (t_pwm_mode *)
        &helicoptercomp_DW.HILInitialize_POModeValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicoptercomp_M, _rt_error_message);
        return;
      }

      {
        int_T i1;
        const uint32_T *p_HILInitialize_POChannels =
          helicoptercomp_P.HILInitialize_POChannels;
        int32_T *dw_POModeValues =
          &helicoptercomp_DW.HILInitialize_POModeValues[0];
        for (i1=0; i1 < 8; i1++) {
          if (dw_POModeValues[i1] == PWM_DUTY_CYCLE_MODE || dw_POModeValues[i1] ==
              PWM_ONE_SHOT_MODE || dw_POModeValues[i1] == PWM_TIME_MODE ||
              dw_POModeValues[i1] == PWM_RAW_MODE) {
            helicoptercomp_DW.HILInitialize_POSortedChans[num_duty_cycle_modes] =
              (p_HILInitialize_POChannels[i1]);
            helicoptercomp_DW.HILInitialize_POSortedFreqs[num_duty_cycle_modes] =
              helicoptercomp_P.HILInitialize_POFrequency;
            num_duty_cycle_modes++;
          } else {
            helicoptercomp_DW.HILInitialize_POSortedChans[7U -
              num_frequency_modes] = (p_HILInitialize_POChannels[i1]);
            helicoptercomp_DW.HILInitialize_POSortedFreqs[7U -
              num_frequency_modes] = helicoptercomp_P.HILInitialize_POFrequency;
            num_frequency_modes++;
          }
        }
      }

      if (num_duty_cycle_modes > 0) {
        result = hil_set_pwm_frequency(helicoptercomp_DW.HILInitialize_Card,
          &helicoptercomp_DW.HILInitialize_POSortedChans[0],
          num_duty_cycle_modes, &helicoptercomp_DW.HILInitialize_POSortedFreqs[0]);
        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(helicoptercomp_M, _rt_error_message);
          return;
        }
      }

      if (num_frequency_modes > 0) {
        result = hil_set_pwm_duty_cycle(helicoptercomp_DW.HILInitialize_Card,
          &helicoptercomp_DW.HILInitialize_POSortedChans[num_duty_cycle_modes],
          num_frequency_modes,
          &helicoptercomp_DW.HILInitialize_POSortedFreqs[num_duty_cycle_modes]);
        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(helicoptercomp_M, _rt_error_message);
          return;
        }
      }

      {
        int_T i1;
        int32_T *dw_POModeValues =
          &helicoptercomp_DW.HILInitialize_POModeValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POModeValues[i1] = helicoptercomp_P.HILInitialize_POConfiguration;
        }
      }

      {
        int_T i1;
        int32_T *dw_POAlignValues =
          &helicoptercomp_DW.HILInitialize_POAlignValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POAlignValues[i1] = helicoptercomp_P.HILInitialize_POAlignment;
        }
      }

      {
        int_T i1;
        int32_T *dw_POPolarityVals =
          &helicoptercomp_DW.HILInitialize_POPolarityVals[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POPolarityVals[i1] = helicoptercomp_P.HILInitialize_POPolarity;
        }
      }

      result = hil_set_pwm_configuration(helicoptercomp_DW.HILInitialize_Card,
        helicoptercomp_P.HILInitialize_POChannels, 8U,
        (t_pwm_configuration *) &helicoptercomp_DW.HILInitialize_POModeValues[0],
        (t_pwm_alignment *) &helicoptercomp_DW.HILInitialize_POAlignValues[0],
        (t_pwm_polarity *) &helicoptercomp_DW.HILInitialize_POPolarityVals[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicoptercomp_M, _rt_error_message);
        return;
      }

      {
        int_T i1;
        real_T *dw_POSortedFreqs =
          &helicoptercomp_DW.HILInitialize_POSortedFreqs[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POSortedFreqs[i1] = helicoptercomp_P.HILInitialize_POLeading;
        }
      }

      {
        int_T i1;
        real_T *dw_POValues = &helicoptercomp_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = helicoptercomp_P.HILInitialize_POTrailing;
        }
      }

      result = hil_set_pwm_deadband(helicoptercomp_DW.HILInitialize_Card,
        helicoptercomp_P.HILInitialize_POChannels, 8U,
        &helicoptercomp_DW.HILInitialize_POSortedFreqs[0],
        &helicoptercomp_DW.HILInitialize_POValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicoptercomp_M, _rt_error_message);
        return;
      }
    }

    if ((helicoptercomp_P.HILInitialize_POStart && !is_switching) ||
        (helicoptercomp_P.HILInitialize_POEnter && is_switching)) {
      {
        int_T i1;
        real_T *dw_POValues = &helicoptercomp_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = helicoptercomp_P.HILInitialize_POInitial;
        }
      }

      result = hil_write_pwm(helicoptercomp_DW.HILInitialize_Card,
        helicoptercomp_P.HILInitialize_POChannels, 8U,
        &helicoptercomp_DW.HILInitialize_POValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicoptercomp_M, _rt_error_message);
        return;
      }
    }

    if (helicoptercomp_P.HILInitialize_POReset) {
      {
        int_T i1;
        real_T *dw_POValues = &helicoptercomp_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = helicoptercomp_P.HILInitialize_POWatchdog;
        }
      }

      result = hil_watchdog_set_pwm_expiration_state
        (helicoptercomp_DW.HILInitialize_Card,
         helicoptercomp_P.HILInitialize_POChannels, 8U,
         &helicoptercomp_DW.HILInitialize_POValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicoptercomp_M, _rt_error_message);
        return;
      }
    }
  }

  /* Start for S-Function (hil_read_encoder_timebase_block): '<S4>/HIL Read Encoder Timebase' */

  /* S-Function Block: helicoptercomp/Helicopter_interface/HIL Read Encoder Timebase (hil_read_encoder_timebase_block) */
  {
    t_error result;
    result = hil_task_create_encoder_reader(helicoptercomp_DW.HILInitialize_Card,
      helicoptercomp_P.HILReadEncoderTimebase_SamplesI,
      helicoptercomp_P.HILReadEncoderTimebase_Channels, 3,
      &helicoptercomp_DW.HILReadEncoderTimebase_Task);
    if (result >= 0) {
      result = hil_task_set_buffer_overflow_mode
        (helicoptercomp_DW.HILReadEncoderTimebase_Task, (t_buffer_overflow_mode)
         (helicoptercomp_P.HILReadEncoderTimebase_Overflow - 1));
    }

    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helicoptercomp_M, _rt_error_message);
    }
  }

  /* Start for FromWorkspace: '<Root>/From Workspace' */
  {
    static real_T pTimeValues0[] = { 0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75,
      2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.25, 4.5, 4.75, 5.0,
      5.25, 5.5, 5.75, 6.0, 6.25, 6.5, 6.75, 7.0, 7.25, 7.5, 7.75, 8.0, 8.25,
      8.5, 8.75, 9.0, 9.25, 9.5, 9.75, 10.0, 10.25, 10.5, 10.75, 11.0, 11.25,
      11.5, 11.75, 12.0, 12.25, 12.5, 12.75, 13.0, 13.25, 13.5, 13.75, 14.0,
      14.25, 14.5, 14.75, 15.0, 15.25, 15.5, 15.75, 16.0, 16.25, 16.5, 16.75,
      17.0, 17.25, 17.5, 17.75, 18.0, 18.25, 18.5, 18.75, 19.0, 19.25, 19.5,
      19.75, 20.0, 20.25, 20.5, 20.75, 21.0, 21.25, 21.5, 21.75, 22.0, 22.25,
      22.5, 22.75, 23.0, 23.25, 23.5, 23.75, 24.0, 24.25, 24.5, 24.75, 25.0,
      25.25, 25.5, 25.75, 26.0, 26.25, 26.5, 26.75, 27.0, 27.25, 27.5, 27.75,
      28.0, 28.25, 28.5, 28.75, 29.0, 29.25, 29.5, 29.75, 30.0, 30.25, 30.5,
      30.75, 31.0, 31.25, 31.5, 31.75, 32.0, 32.25, 32.5, 32.75, 33.0, 33.25,
      33.5, 33.75, 34.0, 34.25, 34.5, 34.75, 35.0 } ;

    static real_T pDataValues0[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.52359877560561874,
      0.52359877560243773, 0.52359877559967882, 0.52359877560176293,
      0.52359877560250567, 0.52359877559879708, 0.523598775601886,
      0.48789589211690976, 0.34971187068198772, 0.22862453939368321,
      0.12361824412228262, 0.033621435678587663, -0.042474371213739148,
      -0.10579616716527829, -0.1574723836113574, -0.19861769233059073,
      -0.23031990327360963, -0.25362901646162456, -0.2695483693600485,
      -0.27902775197638452, -0.28295832332992721, -0.282169144438834,
      -0.2774251374368274, -0.26942628311931704, -0.25880787701330971,
      -0.24614167493421862, -0.2319377716052306, -0.21664706939463596,
      -0.20066420798890355, -0.18433083946977491, -0.16793914653320974,
      -0.1517355142929353, -0.13592427812631835, -0.12067148126441607,
      -0.10610858625019387, -0.092336093964110078, -0.07942703263727624,
      -0.0674302871475923, -0.056373745943501641, -0.046267249192563642,
      -0.0371053272442784, -0.028869723270250169, -0.02153169804541244,
      -0.015054118309729803, -0.0093933330500299839, -0.0045008444163011641,
      -0.00032478188550688625, 0.0031888102427672926, 0.0060948608437053631,
      0.0084483235272315715, 0.0103033793675813, 0.011712778199718943,
      0.012727305885465223, 0.013395365138720638, 0.013762657834992087,
      0.013871957186847195, 0.013762958718040452, 0.013472199591795242,
      0.01303303652261667, 0.012475673208138915, 0.011827228942350554,
      0.011111840800780426, 0.010350792510642193, 0.0095626638251978813,
      0.0087634949041146548, 0.0079669608543352366, 0.007184552204344663,
      0.0064257576653433679, 0.00569824607345959, 0.0050080449064697774,
      0.0043597132260854154, 0.00375650731295063, 0.0032005376369405036,
      0.0026929161415080172, 0.0022338931194832989, 0.0018229832210291397,
      0.0014590803648616424, 0.0011405615241664169, 0.00086537953202614482,
      0.000631145201325567, 0.00043519918545331659, 0.0002746741243856432,
      0.00014654773359501494, 4.7687611652746931E-5, -2.51113174962736E-5,
      -7.5095625605481331E-5, -0.00010552709432021246, -0.00011965954800363665,
      -0.00012071641349342688, -0.00011186568160814758, -9.6187572839689928E-5,
      -7.6628292929425967E-5, -5.5930824064076567E-5, -3.653115928026196E-5,
      -2.0407226918428556E-5, -8.87193000354003E-6, -2.3207706419947627E-6,
      2.2423174428354287E-12, 2.2422064205329661E-12, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0 } ;

    helicoptercomp_DW.FromWorkspace_PWORK.TimePtr = (void *) pTimeValues0;
    helicoptercomp_DW.FromWorkspace_PWORK.DataPtr = (void *) pDataValues0;
    helicoptercomp_DW.FromWorkspace_IWORK.PrevIndex = 0;
  }

  /* Start for If: '<S3>/If' */
  helicoptercomp_DW.If_ActiveSubsystem = -1;

  /* InitializeConditions for TransferFcn: '<S4>/Travel: Transfer Fcn' */
  helicoptercomp_X.TravelTransferFcn_CSTATE = 0.0;

  /* InitializeConditions for TransferFcn: '<S4>/Pitch: Transfer Fcn' */
  helicoptercomp_X.PitchTransferFcn_CSTATE = 0.0;

  /* InitializeConditions for TransferFcn: '<S4>/Elevation: Transfer Fcn' */
  helicoptercomp_X.ElevationTransferFcn_CSTATE = 0.0;

  /* InitializeConditions for Integrator: '<S3>/Integrator' */
  helicoptercomp_X.Integrator_CSTATE = helicoptercomp_P.Integrator_IC;
}

/* Model terminate function */
void helicoptercomp_terminate(void)
{
  /* Terminate for S-Function (hil_initialize_block): '<Root>/HIL Initialize' */

  /* S-Function Block: helicoptercomp/HIL Initialize (hil_initialize_block) */
  {
    t_boolean is_switching;
    t_int result;
    t_uint32 num_final_analog_outputs = 0;
    t_uint32 num_final_pwm_outputs = 0;
    hil_task_stop_all(helicoptercomp_DW.HILInitialize_Card);
    hil_monitor_stop_all(helicoptercomp_DW.HILInitialize_Card);
    is_switching = false;
    if ((helicoptercomp_P.HILInitialize_AOTerminate && !is_switching) ||
        (helicoptercomp_P.HILInitialize_AOExit && is_switching)) {
      {
        int_T i1;
        real_T *dw_AOVoltages = &helicoptercomp_DW.HILInitialize_AOVoltages[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOVoltages[i1] = helicoptercomp_P.HILInitialize_AOFinal;
        }
      }

      num_final_analog_outputs = 8U;
    } else {
      num_final_analog_outputs = 0;
    }

    if ((helicoptercomp_P.HILInitialize_POTerminate && !is_switching) ||
        (helicoptercomp_P.HILInitialize_POExit && is_switching)) {
      {
        int_T i1;
        real_T *dw_POValues = &helicoptercomp_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = helicoptercomp_P.HILInitialize_POFinal;
        }
      }

      num_final_pwm_outputs = 8U;
    } else {
      num_final_pwm_outputs = 0;
    }

    if (0
        || num_final_analog_outputs > 0
        || num_final_pwm_outputs > 0
        ) {
      /* Attempt to write the final outputs atomically (due to firmware issue in old Q2-USB). Otherwise write channels individually */
      result = hil_write(helicoptercomp_DW.HILInitialize_Card
                         , helicoptercomp_P.HILInitialize_AOChannels,
                         num_final_analog_outputs
                         , helicoptercomp_P.HILInitialize_POChannels,
                         num_final_pwm_outputs
                         , NULL, 0
                         , NULL, 0
                         , &helicoptercomp_DW.HILInitialize_AOVoltages[0]
                         , &helicoptercomp_DW.HILInitialize_POValues[0]
                         , (t_boolean *) NULL
                         , NULL
                         );
      if (result == -QERR_HIL_WRITE_NOT_SUPPORTED) {
        t_error local_result;
        result = 0;

        /* The hil_write operation is not supported by this card. Write final outputs for each channel type */
        if (num_final_analog_outputs > 0) {
          local_result = hil_write_analog(helicoptercomp_DW.HILInitialize_Card,
            helicoptercomp_P.HILInitialize_AOChannels, num_final_analog_outputs,
            &helicoptercomp_DW.HILInitialize_AOVoltages[0]);
          if (local_result < 0) {
            result = local_result;
          }
        }

        if (num_final_pwm_outputs > 0) {
          local_result = hil_write_pwm(helicoptercomp_DW.HILInitialize_Card,
            helicoptercomp_P.HILInitialize_POChannels, num_final_pwm_outputs,
            &helicoptercomp_DW.HILInitialize_POValues[0]);
          if (local_result < 0) {
            result = local_result;
          }
        }

        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(helicoptercomp_M, _rt_error_message);
        }
      }
    }

    hil_task_delete_all(helicoptercomp_DW.HILInitialize_Card);
    hil_monitor_delete_all(helicoptercomp_DW.HILInitialize_Card);
    hil_close(helicoptercomp_DW.HILInitialize_Card);
    helicoptercomp_DW.HILInitialize_Card = NULL;
  }
}

/*========================================================================*
 * Start of Classic call interface                                        *
 *========================================================================*/

/* Solver interface called by GRT_Main */
#ifndef USE_GENERATED_SOLVER

void rt_ODECreateIntegrationData(RTWSolverInfo *si)
{
  UNUSED_PARAMETER(si);
  return;
}                                      /* do nothing */

void rt_ODEDestroyIntegrationData(RTWSolverInfo *si)
{
  UNUSED_PARAMETER(si);
  return;
}                                      /* do nothing */

void rt_ODEUpdateContinuousStates(RTWSolverInfo *si)
{
  UNUSED_PARAMETER(si);
  return;
}                                      /* do nothing */

#endif

void MdlOutputs(int_T tid)
{
  if (tid == 1)
    tid = 0;
  helicoptercomp_output(tid);
}

void MdlUpdate(int_T tid)
{
  if (tid == 1)
    tid = 0;
  helicoptercomp_update(tid);
}

void MdlInitializeSizes(void)
{
}

void MdlInitializeSampleTimes(void)
{
}

void MdlInitialize(void)
{
}

void MdlStart(void)
{
  helicoptercomp_initialize();
}

void MdlTerminate(void)
{
  helicoptercomp_terminate();
}

/* Registration function */
RT_MODEL_helicoptercomp_T *helicoptercomp(void)
{
  /* Registration code */

  /* initialize non-finites */
  rt_InitInfAndNaN(sizeof(real_T));

  /* non-finite (run-time) assignments */
  helicoptercomp_P.Integrator_UpperSat = rtInf;
  helicoptercomp_P.Integrator_LowerSat = rtMinusInf;

  /* initialize real-time model */
  (void) memset((void *)helicoptercomp_M, 0,
                sizeof(RT_MODEL_helicoptercomp_T));

  {
    /* Setup solver object */
    rtsiSetSimTimeStepPtr(&helicoptercomp_M->solverInfo,
                          &helicoptercomp_M->Timing.simTimeStep);
    rtsiSetTPtr(&helicoptercomp_M->solverInfo, &rtmGetTPtr(helicoptercomp_M));
    rtsiSetStepSizePtr(&helicoptercomp_M->solverInfo,
                       &helicoptercomp_M->Timing.stepSize0);
    rtsiSetdXPtr(&helicoptercomp_M->solverInfo, &helicoptercomp_M->derivs);
    rtsiSetContStatesPtr(&helicoptercomp_M->solverInfo, (real_T **)
                         &helicoptercomp_M->contStates);
    rtsiSetNumContStatesPtr(&helicoptercomp_M->solverInfo,
      &helicoptercomp_M->Sizes.numContStates);
    rtsiSetNumPeriodicContStatesPtr(&helicoptercomp_M->solverInfo,
      &helicoptercomp_M->Sizes.numPeriodicContStates);
    rtsiSetPeriodicContStateIndicesPtr(&helicoptercomp_M->solverInfo,
      &helicoptercomp_M->periodicContStateIndices);
    rtsiSetPeriodicContStateRangesPtr(&helicoptercomp_M->solverInfo,
      &helicoptercomp_M->periodicContStateRanges);
    rtsiSetErrorStatusPtr(&helicoptercomp_M->solverInfo, (&rtmGetErrorStatus
      (helicoptercomp_M)));
    rtsiSetRTModelPtr(&helicoptercomp_M->solverInfo, helicoptercomp_M);
  }

  rtsiSetSimTimeStep(&helicoptercomp_M->solverInfo, MAJOR_TIME_STEP);
  helicoptercomp_M->intgData.f[0] = helicoptercomp_M->odeF[0];
  helicoptercomp_M->contStates = ((real_T *) &helicoptercomp_X);
  rtsiSetSolverData(&helicoptercomp_M->solverInfo, (void *)
                    &helicoptercomp_M->intgData);
  rtsiSetSolverName(&helicoptercomp_M->solverInfo,"ode1");

  /* Initialize timing info */
  {
    int_T *mdlTsMap = helicoptercomp_M->Timing.sampleTimeTaskIDArray;
    mdlTsMap[0] = 0;
    mdlTsMap[1] = 1;
    mdlTsMap[2] = 2;
    helicoptercomp_M->Timing.sampleTimeTaskIDPtr = (&mdlTsMap[0]);
    helicoptercomp_M->Timing.sampleTimes =
      (&helicoptercomp_M->Timing.sampleTimesArray[0]);
    helicoptercomp_M->Timing.offsetTimes =
      (&helicoptercomp_M->Timing.offsetTimesArray[0]);

    /* task periods */
    helicoptercomp_M->Timing.sampleTimes[0] = (0.0);
    helicoptercomp_M->Timing.sampleTimes[1] = (0.002);
    helicoptercomp_M->Timing.sampleTimes[2] = (0.25);

    /* task offsets */
    helicoptercomp_M->Timing.offsetTimes[0] = (0.0);
    helicoptercomp_M->Timing.offsetTimes[1] = (0.0);
    helicoptercomp_M->Timing.offsetTimes[2] = (0.0);
  }

  rtmSetTPtr(helicoptercomp_M, &helicoptercomp_M->Timing.tArray[0]);

  {
    int_T *mdlSampleHits = helicoptercomp_M->Timing.sampleHitArray;
    int_T *mdlPerTaskSampleHits =
      helicoptercomp_M->Timing.perTaskSampleHitsArray;
    helicoptercomp_M->Timing.perTaskSampleHits = (&mdlPerTaskSampleHits[0]);
    mdlSampleHits[0] = 1;
    helicoptercomp_M->Timing.sampleHits = (&mdlSampleHits[0]);
  }

  rtmSetTFinal(helicoptercomp_M, 35.0);
  helicoptercomp_M->Timing.stepSize0 = 0.002;
  helicoptercomp_M->Timing.stepSize1 = 0.002;
  helicoptercomp_M->Timing.stepSize2 = 0.25;
  helicoptercomp_M->solverInfoPtr = (&helicoptercomp_M->solverInfo);
  helicoptercomp_M->Timing.stepSize = (0.002);
  rtsiSetFixedStepSize(&helicoptercomp_M->solverInfo, 0.002);
  rtsiSetSolverMode(&helicoptercomp_M->solverInfo, SOLVER_MODE_MULTITASKING);

  /* block I/O */
  helicoptercomp_M->blockIO = ((void *) &helicoptercomp_B);

  {
    helicoptercomp_B.TravelCounttorad = 0.0;
    helicoptercomp_B.Sum3 = 0.0;
    helicoptercomp_B.PitchCounttorad = 0.0;
    helicoptercomp_B.Gain = 0.0;
    helicoptercomp_B.ElevationCounttorad = 0.0;
    helicoptercomp_B.Sum = 0.0;
    helicoptercomp_B.FrontmotorSaturation = 0.0;
    helicoptercomp_B.BackmotorSaturation = 0.0;
    helicoptercomp_B.In1 = 0.0;
  }

  /* parameters */
  helicoptercomp_M->defaultParam = ((real_T *)&helicoptercomp_P);

  /* states (continuous) */
  {
    real_T *x = (real_T *) &helicoptercomp_X;
    helicoptercomp_M->contStates = (x);
    (void) memset((void *)&helicoptercomp_X, 0,
                  sizeof(X_helicoptercomp_T));
  }

  /* states (dwork) */
  helicoptercomp_M->dwork = ((void *) &helicoptercomp_DW);
  (void) memset((void *)&helicoptercomp_DW, 0,
                sizeof(DW_helicoptercomp_T));

  {
    int32_T i;
    for (i = 0; i < 8; i++) {
      helicoptercomp_DW.HILInitialize_AIMinimums[i] = 0.0;
    }
  }

  {
    int32_T i;
    for (i = 0; i < 8; i++) {
      helicoptercomp_DW.HILInitialize_AIMaximums[i] = 0.0;
    }
  }

  {
    int32_T i;
    for (i = 0; i < 8; i++) {
      helicoptercomp_DW.HILInitialize_AOMinimums[i] = 0.0;
    }
  }

  {
    int32_T i;
    for (i = 0; i < 8; i++) {
      helicoptercomp_DW.HILInitialize_AOMaximums[i] = 0.0;
    }
  }

  {
    int32_T i;
    for (i = 0; i < 8; i++) {
      helicoptercomp_DW.HILInitialize_AOVoltages[i] = 0.0;
    }
  }

  {
    int32_T i;
    for (i = 0; i < 8; i++) {
      helicoptercomp_DW.HILInitialize_FilterFrequency[i] = 0.0;
    }
  }

  {
    int32_T i;
    for (i = 0; i < 8; i++) {
      helicoptercomp_DW.HILInitialize_POSortedFreqs[i] = 0.0;
    }
  }

  {
    int32_T i;
    for (i = 0; i < 8; i++) {
      helicoptercomp_DW.HILInitialize_POValues[i] = 0.0;
    }
  }

  helicoptercomp_DW.HILWriteAnalog_Buffer[0] = 0.0;
  helicoptercomp_DW.HILWriteAnalog_Buffer[1] = 0.0;

  /* Initialize Sizes */
  helicoptercomp_M->Sizes.numContStates = (4);/* Number of continuous states */
  helicoptercomp_M->Sizes.numPeriodicContStates = (0);
                                      /* Number of periodic continuous states */
  helicoptercomp_M->Sizes.numY = (0);  /* Number of model outputs */
  helicoptercomp_M->Sizes.numU = (0);  /* Number of model inputs */
  helicoptercomp_M->Sizes.sysDirFeedThru = (0);/* The model is not direct feedthrough */
  helicoptercomp_M->Sizes.numSampTimes = (3);/* Number of sample times */
  helicoptercomp_M->Sizes.numBlocks = (48);/* Number of blocks */
  helicoptercomp_M->Sizes.numBlockIO = (9);/* Number of block outputs */
  helicoptercomp_M->Sizes.numBlockPrms = (144);/* Sum of parameter "widths" */
  return helicoptercomp_M;
}

/*========================================================================*
 * End of Classic call interface                                          *
 *========================================================================*/
