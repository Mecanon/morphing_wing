
const int transistorPin = 9;    // connected to the base of the transistor

 void setup() {
   // set  the transistor pin as output:
   pinMode(transistorPin, OUTPUT);
   Serial.begin(9600);
 }
 
 void loop() {
  for(int i = 0; i<360; i++){
    //convert 0-360 angle to radian (needed for sin function)
    float rad = DEG_TO_RAD * i;
   //calculate sin of angle as number between 0 and 255
   int sinOut = constrain((sin(rad) * 128) + 128, 0, 255); 
   // map the sensor value to a range from 0 - 255:
   // int outputValue = map(sensorValue, 0, 1023, 0, 255);
   // use that to control the transistor:
   analogWrite(transistorPin, sinOut);
   Serial.println(sinOut);
   delay(50);
  }
 }
