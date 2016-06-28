
/*
Original Author:
Basic Example Sketch for the ITG-3200 (http://www.sparkfun.com/products/9801)
SparkFun Electronics 2011
Ryan Owens
3/16/11
Updtaed for Arduino 1.0 and beyond:
Joel Bartlett
SparkFun Electronics
10/16/13
This code is public domain buy you buy me a beer if you use this and we meet someday (Beerware license).
To use this example code, attach:
Arduino  :  ITG-3200 Breakout
3.3V  :  VDD
3.3V  :  VIO
GND   :  GND
SDA   :  SDA or A4 on older Arduinos
SCL   :  SCL or A5 on older Arduinos
Load the sketch and open the serial window at 9600 bps. Arduino will output the raw X,Y and Z axis data being read from the gyroscope.
*/

//The Wire library is used for I2C communication
#include <Wire.h>
#include <RunningMedian.h>


RunningMedian calibration_data_x = RunningMedian(100);
RunningMedian calibration_data_y = RunningMedian(100);
RunningMedian calibration_data_z = RunningMedian(100);

//This is a list of registers in the ITG-3200. Registers are parameters that determine how the sensor will behave, or they can hold data that represent the
//sensors current status.
//To learn more about the registers on the ITG-3200, download and read the datasheet.
char WHO_AM_I = 0x00;
char SMPLRT_DIV = 0x15;
char DLPF_FS = 0x16;
char GYRO_XOUT_H = 0x1D;
char GYRO_XOUT_L = 0x1E;
char GYRO_YOUT_H = 0x1F;
char GYRO_YOUT_L = 0x20;
char GYRO_ZOUT_H = 0x21;
char GYRO_ZOUT_L = 0x22;

//This is a list of settings that can be loaded into the registers.
//DLPF, Full Scale Register Bits
//FS_SEL must be set to 3 for proper operation
//Set DLPF_CFG to 3 for 1kHz Fint and 42 Hz Low Pass Filter
char DLPF_CFG_0 = (1 << 0);
char DLPF_CFG_1 = (1 << 1);
char DLPF_CFG_2 = (1 << 2);
char DLPF_FS_SEL_0 = (1 << 3);
char DLPF_FS_SEL_1 = (1 << 4);

//I2C devices each have an address. The address is defined in the datasheet for the device. The ITG-3200 breakout board can have different address depending on how
//the jumper on top of the board is configured. By default, the jumper is connected to the VDD pin. When the jumper is connected to the VDD pin the I2C address
//is 0x69.
char itgAddress = 0x69;

//Values to be callibrated
float x_cal, y_cal, z_cal;
float x_angle = 0.;
float y_angle = 0.;
float z_angle = 0.;
float calibration = 14.375;
int step_time = 50;

const int transistorPin = 9;    // connected to the base of the transistor
int i = 0;

//In the setup section of the sketch the serial port will be configured, the i2c communication will be initialized, and the itg-3200 will be configured.
void setup()
{
  // set  the transistor pin as output:
  pinMode(transistorPin, OUTPUT);
  //Create a serial connection using a 9600bps baud rate.
  Serial.begin(9600);

  //Initialize the I2C communication. This will set the Arduino up as the 'Master' device.
  Wire.begin();

  //Read the WHO_AM_I register and print the result
  char id = 0;
  id = itgRead(itgAddress, 0x00);
  Serial.print("ID: ");
  Serial.println(id, HEX);

  //Configure the gyroscope
  //Set the gyroscope scale for the outputs to +/-2000 degrees per second
  itgWrite(itgAddress, DLPF_FS, (DLPF_FS_SEL_0 | DLPF_FS_SEL_1 | DLPF_CFG_0));
  //Set the sample rate to 100 hz
  itgWrite(itgAddress, SMPLRT_DIV, 9);

  //Calibrate values

  for (int j = 0; j < 20; j++) {
    calibration_data_x.add(readX());
    calibration_data_y.add(readY());
    calibration_data_z.add(readZ());
    delay(step_time);
  }
  x_cal = calibration_data_x.getAverage() / calibration;
  y_cal = calibration_data_y.getAverage() / calibration;
  z_cal = calibration_data_z.getAverage() / calibration;
  //  Serial.print("Calibrated values:");
  //  Serial.print(x_cal);
  //  Serial.print('\t');
  //  Serial.print(y_cal);
  //  Serial.print('\t');
  //  Serial.println(z_cal);
  //Serial.print('\n');
}

//The loop section of the sketch will read the X,Y and Z output rates from the gyroscope and output them in the Serial Terminal
void loop()
{
  //convert 0-360 angle to radian (needed for sin function)
  float rad = DEG_TO_RAD * i;
  //calculate sin of angle as number between 0 and 255
  int sinOut = constrain((sin(rad) * 128) + 128, 0, 255);
  // map the sensor value to a range from 0 - 255:
  // int outputValue = map(sensorValue, 0, 1023, 0, 255);
  // use that to control the transistor:
  analogWrite(transistorPin, sinOut);
  Serial.print(sinOut/255.);
  Serial.print('\t');
  i = i + 1;

  //Create variables to hold the output rates.
  float xRate, yRate, zRate;

  //Read the x,y and z output rates from the gyroscope.
  xRate = readX() / calibration - x_cal;
  yRate = readY() / calibration - y_cal;
  zRate = readZ() / calibration - z_cal;

  x_angle = x_angle + xRate * step_time / 1000.;
  y_angle = y_angle + yRate * step_time / 1000.;
  z_angle = z_angle + zRate * step_time / 1000.;

  //Print the output rates to the terminal, seperated by a TAB character.
  Serial.print(x_angle);
  Serial.print('\t');
  Serial.print(y_angle);
  Serial.print('\t');
  Serial.println(z_angle);
  // Serial.print('\n');
  //Wait 10ms before reading the values again. (Remember, the output rate was set to 100hz and 1reading per 10ms = 100hz.)
  delay(step_time);
}

//This function will write a value to a register on the itg-3200.
//Parameters:
//  char address: The I2C address of the sensor. For the ITG-3200 breakout the address is 0x69.
//  char registerAddress: The address of the register on the sensor that should be written to.
//  char data: The value to be written to the specified register.
void itgWrite(char address, char registerAddress, char data)
{
  //Initiate a communication sequence with the desired i2c device
  Wire.beginTransmission(address);
  //Tell the I2C address which register we are writing to
  Wire.write(registerAddress);
  //Send the value to write to the specified register
  Wire.write(data);
  //End the communication sequence
  Wire.endTransmission();
}

//This function will read the data from a specified register on the ITG-3200 and return the value.
//Parameters:
//  char address: The I2C address of the sensor. For the ITG-3200 breakout the address is 0x69.
//  char registerAddress: The address of the register on the sensor that should be read
//Return:
//  unsigned char: The value currently residing in the specified register
unsigned char itgRead(char address, char registerAddress)
{
  //This variable will hold the contents read from the i2c device.
  unsigned char data = 0;

  //Send the register address to be read.
  Wire.beginTransmission(address);
  //Send the Register Address
  Wire.write(registerAddress);
  //End the communication sequence.
  Wire.endTransmission();

  //Ask the I2C device for data
  Wire.beginTransmission(address);
  Wire.requestFrom(address, 1);

  //Wait for a response from the I2C device
  if (Wire.available()) {
    //Save the data sent from the I2C device
    data = Wire.read();
  }

  //End the communication sequence.
  Wire.endTransmission();

  //Return the data read during the operation
  return data;
}

//This function is used to read the X-Axis rate of the gyroscope. The function returns the ADC value from the Gyroscope
//NOTE: This value is NOT in degrees per second.
//Usage: int xRate = readX();
int readX(void)
{
  int data = 0;
  data = itgRead(itgAddress, GYRO_XOUT_H) << 8;
  data |= itgRead(itgAddress, GYRO_XOUT_L);

  return data;
}

//This function is used to read the Y-Axis rate of the gyroscope. The function returns the ADC value from the Gyroscope
//NOTE: This value is NOT in degrees per second.
//Usage: int yRate = readY();
int readY(void)
{
  int data = 0;
  data = itgRead(itgAddress, GYRO_YOUT_H) << 8;
  data |= itgRead(itgAddress, GYRO_YOUT_L);

  return data;
}

//This function is used to read the Z-Axis rate of the gyroscope. The function returns the ADC value from the Gyroscope
//NOTE: This value is NOT in degrees per second.
//Usage: int zRate = readZ();
int readZ(void)
{
  int data = 0;
  data = itgRead(itgAddress, GYRO_ZOUT_H) << 8;
  data |= itgRead(itgAddress, GYRO_ZOUT_L);

  return data;
}
