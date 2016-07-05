int analogInPin = A0;  // Analog input pin that the carrier board OUT is connected to

double sensorValue = 0;        // value read from the carrier board
double outputValue = 0;        // output in milliamps
int mVperAmp = 185*4; // use 100 for 20A Module and 66 for 30A Module
int ACSoffset = 2500; //2500; 
double Voltage = 0;
double instCurrent =0;
double Current =0;
double Value =0;
double Media = 0;
const int transistorPin = 9;    // connected to the base of the transistor

double inputVoltage = 3000; // Put here the voltage shown in the sorce 

void setup() {
  // initialize serial communications at 9600 bps:
  Serial.begin(9600);
  pinMode(transistorPin, OUTPUT);

}

void loop() {
  instCurrent = 0;
  for(int i = 0; i<3600; i++){
    //convert 0-360 angle to radian (needed for sin function)
    float rad = DEG_TO_RAD * i/10;
   //calculate sin of angle as number between 0 and 255
   int sinOut = 255; //constrain((cos(rad - 3.1415) * 128) + 128, 0, 255); 
   // map the sensor value to a range from 0 - 255:
    //int outputValue = map(sensorValue, 0, 1023, 0, 255);
   // use that to control the transistor:
   analogWrite(transistorPin, sinOut);
   Value = (inputVoltage*sinOut)/255.0;
   Serial.print(Value);
   Serial.print("  ");
   //delay(50);
   
  // read the analog in value:
  sensorValue = analogRead(analogInPin);
  // convert to milli amps
  Current = ((((sensorValue / 1024.0) * 5000) - ACSoffset) / mVperAmp)*1000;
  Serial.print(Current);
  Serial.println(" ");
  

  // wait 10 milliseconds before the next loop
  // for the analog-to-digital converter to settle
  // after the last reading:
  delay(1000);
  }

}
